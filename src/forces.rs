use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use crate::particle::{GravitationalSource, Particle};
use crate::gravitree::{Node, AccuracyCriterion};
use crate::vectors::*;

pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const KM_IN_KPC: f64 = 30856776000000000.0; // number of km in kpc

#[derive(Debug)]
#[derive(Clone)]
pub enum ForceCalculationMethod {
    Direct,
    Tree(AccuracyCriterion),
}

impl ForceCalculationMethod {
    pub fn name(&self) -> String {
        match self {
            ForceCalculationMethod::Direct => String::from("Direct"),
            ForceCalculationMethod::Tree(criterion) => String::from("Tree") + "-" + &criterion.name()
        }
    }
}

pub fn recalculate_dynamics_due_to_gravity_directly(data: &mut Vec<Particle>) {
    for i in 0..data.len() {
        // Save for use with jerk and dynamical tree criterion
        let previous_acceleration = data[i].acceleration.clone();
        
        // reset acceleration, jerk, snap
        data[i].acceleration = [0.0; 3];
        data[i].jerk = [0.0; 3];
        data[i].snap = [0.0; 3];

        // calculate dynamics via direct summation
        for j in 0..data.len() {
            if i == j {
                continue;
            };
            let source = data[j].clone();
            calculate_dynamics_due_to_one_body(&mut data[i], &previous_acceleration, &source);
        }
    }
}

pub fn recalculate_dynamics_due_to_gravity_with_tree(data: &mut Vec<Particle>, criterion: &AccuracyCriterion, root: &Box<Node>) {
    data.par_iter_mut().for_each(|particle| {
        // Save for use with jerk and dynamical tree criterion
        let previous_acceleration = particle.acceleration.clone();
        
        // reset acceleration, jerk, snap
        particle.acceleration = [0.0; 3];
        particle.jerk = [0.0; 3];
        particle.snap = [0.0; 3];
        
        // calculate dynamics via tree
        calculate_dynamics_with_tree(particle, &previous_acceleration, root, &criterion);
    });
}

pub fn calculate_dynamics_with_tree(target: &mut Particle, previous_target_accel: &[f64; 3], node: &Box<Node>, criterion: &AccuracyCriterion) {
    match &**node {
        Node::Leaf { particle } => {
            if node.contains(target) {
                return;
            }
        
            calculate_dynamics_due_to_one_body(target, &previous_target_accel, particle);
        },
        Node::Branch { children, .. } => {
            if node.is_approximatable(target, previous_target_accel, criterion) && !node.contains(&target) {
                calculate_dynamics_due_to_one_body(target, &previous_target_accel, &**node);
            } else {
                for child in children.iter() {
                    calculate_dynamics_with_tree(target, previous_target_accel, child, criterion);
                }
            }
        }
    }
}

fn calculate_dynamics_due_to_one_body<T: GravitationalSource + std::fmt::Debug>(target: &mut Particle, previous_target_accel: &[f64; 3], source: &T) {
    let displacement: [f64; 3] = subtract(&source.get_position(), &target.position);
    let displacement_mag: f64 = magnitude(&displacement);
    let inv_r3 = displacement_mag.powi(-3);
    let inv_r2 = displacement_mag.powi(-2);
    let inv_r5 = inv_r2 * inv_r3;
    let inv_r7 = inv_r5 * inv_r2;
    if displacement_mag == 0.0 {
        panic!("Zero displacement in dynamics calculation! This probably means a particle is errorously interacting with itself.");
    }

    let relative_velocity: [f64; 3] = subtract(&source.get_velocity(), &target.velocity);
    let relative_velocity_squared: f64 = dot_product(&relative_velocity, &relative_velocity);
    let rdotv: f64 = dot_product(&displacement, &relative_velocity);

    let mut relative_accel: [f64; 3] = subtract(&source.get_acceleration(), &previous_target_accel); // puts in units of km^2 / kpc s^2
    let rdota: f64 = dot_product(&displacement, &relative_accel) * KM_IN_KPC; // puts it in units of km^2 / s^2

    let mut acceleration: [f64; 3] = [0.0; 3];
    let mut jerk: [f64; 3] = [0.0; 3];
    let mut snap: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        // put relative accel in km^2 / kpc s^2 so it has same units as other terms in snap sum
        relative_accel[i] *= KM_IN_KPC;

        // units of accel are naively km^2 / kpc s^2 but we want km / s^2
        acceleration[i] = GG * source.get_mass() * displacement[i] * inv_r3;
        acceleration[i] /= KM_IN_KPC;
    
        // units of jerk are naively km^3 / kpc^2 s^2 but we want km / s^3
        jerk[i] = GG * source.get_mass() * ((relative_velocity[i] * inv_r3) - 3.0 * rdotv * (displacement[i] * inv_r5));
        jerk[i] /= KM_IN_KPC.powi(2);

        // units of snap are naively km^4 / kpc^3 s^3 but we want km / s^4
        snap[i] = GG * source.get_mass() *(
            (relative_accel[i] * inv_r3)
            - 6.0 * (rdotv * relative_velocity[i] * inv_r5)
            - 3.0 * (relative_velocity_squared + rdota) * (displacement[i] * inv_r5)
            + 15.0 * rdotv.powi(2) * (displacement[i] * inv_r7)
        );
        snap[i] /= KM_IN_KPC.powi(3);
    };

    // Apply changes
    target.acceleration = add(&target.acceleration, &acceleration);
    target.jerk = add(&target.jerk, &jerk);
    target.snap = add(&target.snap, &snap);
}

pub fn calculate_energy(data: &Vec<Particle>) -> (f64, f64) {
    // units of energy are Msun km^2 / s^2
    let mut potential: f64 = 0.0;
    let mut kinetic: f64 = 0.0;
    for i in 0..data.len() {
        kinetic += kinetic_energy(&data[i]);

        for j in 0..data.len() {
            if i <= j {
                continue;
            };
            potential += potential_of_binary_config(&data[i], &data[j]);
        };
    };
    (kinetic, potential)
}

fn kinetic_energy(particle: &Particle) -> f64 {
    0.5 * particle.mass * dot_product(&particle.velocity, &particle.velocity)
}

fn potential_of_binary_config(particle1: &Particle, particle2: &Particle) -> f64 {
    -1.0 * GG * particle1.mass * particle2.mass / magnitude(&subtract(&particle2.position, &particle1.position))
}

pub fn calculate_momentum(data: &Vec<Particle>) -> [f64; 3] {
    let mut momentum: [f64; 3] = [0.0; 3];
    for i in 0..data.len() {
        momentum = add(&momentum, &scalar_multiply(&data[i].mass, &data[i].velocity))
    }
    momentum
}

pub fn calculate_angular_momentum(data: &Vec<Particle>) -> [f64; 3] {
    // units of energy are Msun km kpc / s
    let mut angular_momentum: [f64; 3] = [0.0, 0.0, 0.0];
    for i in 0..data.len() {
        let momentum = scalar_multiply(&data[i].mass, &data[i].velocity);
        angular_momentum = add(&angular_momentum, &cross_product(&data[i].position, &momentum));
    };
    angular_momentum
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;
    use crate::time_evol::compute_timestep;
    use crate::init_conds::plummer_init_conds;
    use crate::gravitree::build_gravitree;

    #[test]
    fn test_dynamics() {
        let mut earth = Particle {
            mass: 3.0e-6,
            position: [4.848136811e-9, 0.0, 0.0],
            velocity: [0.0, 29.78, 0.0],
            acceleration: [-5.930e-6, 0.0, 0.0],
            jerk: [0.0; 3],
            snap: [0.0; 3]
        };

        let sun = Particle {
            mass: 1.0,
            position: [0.0; 3],
            velocity: [0.0; 3],
            acceleration: [0.0; 3],
            jerk: [0.0; 3],
            snap: [0.0; 3]
        };

        let previous_earth_accel: [f64; 3] = earth.acceleration.clone();

        //reset dynamics
        earth.acceleration = [0.0; 3];
        earth.jerk = [0.0; 3];
        earth.snap = [0.0; 3];

        calculate_dynamics_due_to_one_body(&mut earth, &previous_earth_accel, &sun);

        let timestep = compute_timestep(&earth, &0.1);
        println!("dt = {}", timestep);

        let expected_acceleration = 5.930e-6;
        let expected_jerk = 1.1805e-12;
        let expected_snap = 2.3484e-19;
        let expected_timestep = 3.55260e5;

        let accel_error = (magnitude(&earth.acceleration)/expected_acceleration) - 1.0;
        let jerk_error = (magnitude(&earth.jerk)/expected_jerk) - 1.0;
        let snap_error = (magnitude(&earth.snap)/expected_snap) - 1.0;
        let timestep_error = (timestep/expected_timestep) - 1.0;

        if accel_error.abs() > 1e-4 {
            panic!("Acceleration calculation is inaccurate! Error = {}", accel_error)
        };

        if jerk_error.abs() > 1e-4 {
            panic!("Jerk calculation is inaccurate! Error = {}", jerk_error)
        };

        if snap_error.abs() > 1e-4 {
            panic!("Snap calculation is inaccurate! Error = {}", snap_error)
        };

        if timestep_error.abs() > 1e-4 {
            panic!("Timestep calculation is inaccurate! Error = {}", timestep_error)
        };

    }

    #[test]
    fn test_tree() {
        let mut particles = plummer_init_conds(1.0, 1e8, 1000, String::from("tests"));

        recalculate_dynamics_due_to_gravity_directly(&mut particles);
        let start = Instant::now();
        recalculate_dynamics_due_to_gravity_directly(&mut particles); //two calls to "warm up" higher order derivatives
        let duration = start.elapsed();
        println!("Direct calculation time for {} particles: {:?}", particles.len(), duration);

        let direct_accels: Vec<[f64; 3]> = particles.iter().map(|p| p.acceleration).collect();
        let direct_jerks: Vec<[f64; 3]> = particles.iter().map(|p| p.jerk).collect();
        let direct_snaps: Vec<[f64; 3]> = particles.iter().map(|p| p.snap).collect();

        let start = Instant::now();
        let root = Box::new(build_gravitree(particles.clone()));
        recalculate_dynamics_due_to_gravity_with_tree(&mut particles, &AccuracyCriterion::Geometric(0.3), &root);
        let duration = start.elapsed();
        println!("Geometric tree calculation time for {} particles: {:?}", particles.len(), duration);

        let geometric_tree_accels: Vec<[f64; 3]> = particles.iter().map(|p| p.acceleration).collect();
        let geometric_tree_jerks: Vec<[f64; 3]> = particles.iter().map(|p| p.jerk).collect();
        let geometric_tree_snaps: Vec<[f64; 3]> = particles.iter().map(|p| p.snap).collect();

        let start = Instant::now();
        let root = Box::new(build_gravitree(particles.clone()));
        recalculate_dynamics_due_to_gravity_with_tree(&mut particles, &AccuracyCriterion::Dynamical(1e-3), &root);
        let duration = start.elapsed();
        println!("Dynamical tree calculation time for {} particles: {:?}", particles.len(), duration);

        let dynamical_tree_accels: Vec<[f64; 3]> = particles.iter().map(|p| p.acceleration).collect();
        let dynamical_tree_jerks: Vec<[f64; 3]> = particles.iter().map(|p| p.jerk).collect();
        let dynamical_tree_snaps: Vec<[f64; 3]> = particles.iter().map(|p| p.snap).collect();

        let (geo_accel_err, geo_jerk_err, geo_snap_err) = {
            let mut rms_accel_err = 0.0;
            let mut rms_jerk_err = 0.0;
            let mut rms_snap_err = 0.0;
            for i in 0..particles.len() {
                let accel_err = magnitude(&subtract(&geometric_tree_accels[i], &direct_accels[i])) / magnitude(&direct_accels[i]);
                let jerk_err = magnitude(&subtract(&geometric_tree_jerks[i], &direct_jerks[i])) / magnitude(&direct_jerks[i]);
                let snap_err = magnitude(&subtract(&geometric_tree_snaps[i], &direct_snaps[i])) / magnitude(&direct_snaps[i]);

                rms_accel_err += accel_err * accel_err;
                rms_jerk_err += jerk_err * jerk_err;
                rms_snap_err += snap_err * snap_err;
            }
            (
                rms_accel_err.sqrt() / (particles.len() as f64).sqrt(),
                rms_jerk_err.sqrt() / (particles.len() as f64).sqrt(),
                rms_snap_err.sqrt() / (particles.len() as f64).sqrt()
            )
        };

        let (dyn_accel_err, dyn_jerk_err, dyn_snap_err) = {
            let mut rms_accel_err = 0.0;
            let mut rms_jerk_err = 0.0;
            let mut rms_snap_err = 0.0;
            for i in 0..particles.len() {
                let accel_err = magnitude(&subtract(&dynamical_tree_accels[i], &direct_accels[i])) / magnitude(&direct_accels[i]);
                let jerk_err = magnitude(&subtract(&dynamical_tree_jerks[i], &direct_jerks[i])) / magnitude(&direct_jerks[i]);
                let snap_err = magnitude(&subtract(&dynamical_tree_snaps[i], &direct_snaps[i])) / magnitude(&direct_snaps[i]);

                rms_accel_err += accel_err * accel_err;
                rms_jerk_err += jerk_err * jerk_err;
                rms_snap_err += snap_err * snap_err;
            }
            (
                rms_accel_err.sqrt() / (particles.len() as f64).sqrt(),
                rms_jerk_err.sqrt() / (particles.len() as f64).sqrt(),
                rms_snap_err.sqrt() / (particles.len() as f64).sqrt()
            )
        };

        if geo_accel_err > 1e-2 {
            panic!("Geometric tree acceleration error too high! RMS relative error = {}", geo_accel_err)
        };

        if geo_jerk_err > 5e-2 {
            panic!("Geometric tree jerk error too high! RMS relative error = {}", geo_jerk_err)
        };

        if geo_snap_err > 1e-1 {
            panic!("Geometric tree snap error too high! RMS relative error = {}", geo_snap_err)
        };

        if dyn_accel_err > 1e-2 {
            panic!("Dynamical tree acceleration error too high! RMS relative error = {}", dyn_accel_err)
        };

        if dyn_jerk_err > 5e-2 {
            panic!("Dynamical tree jerk error too high! RMS relative error = {}", dyn_jerk_err)
        };

        if dyn_snap_err > 1e-1 {
            panic!("Dynamical tree snap error too high! RMS relative error = {}", dyn_snap_err)
        };
    }
}