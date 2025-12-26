use crate::particle::Particle;
use crate::vectors::*;

pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const KM_IN_KPC: f64 = 30856776000000000.0; // number of km in kpc

pub fn recalculate_dynamics_due_to_gravity(data: &mut Vec<Particle>) {
    for i in 0..data.len() {
        let previous_acceleration = data[i].acceleration.clone();
        data[i].acceleration = [0.0; 3];
        data[i].jerk = [0.0; 3];
        data[i].snap = [0.0; 3];
        for j in 0..data.len() {
            if i == j {
                continue;
            };
            let dynamical_variables = &calculate_dynamics_due_to_one_body(&data[i], &previous_acceleration, &data[j]);
            data[i].acceleration = add(&data[i].acceleration, &dynamical_variables.0);
            data[i].jerk = add(&data[i].jerk, &dynamical_variables.1);
            data[i].snap = add(&data[i].snap, &dynamical_variables.2);
        };
    };
}

fn calculate_dynamics_due_to_one_body(target: &Particle, previous_target_accel: &[f64; 3], source: &Particle) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let displacement: [f64; 3] = subtract(&source.position, &target.position);
    let displacement_mag: f64 = magnitude(&displacement);
    
    let relative_velocity: [f64; 3] = subtract(&source.velocity, &target.velocity);
    let relative_velocity_squared: f64 = dot_product(&relative_velocity, &relative_velocity);
    let rdotv: f64 = dot_product(&displacement, &relative_velocity);

    let mut relative_accel: [f64; 3] = subtract(&source.acceleration, &previous_target_accel); // puts in units of km^2 / kpc s^2
    let rdota: f64 = dot_product(&displacement, &relative_accel) * KM_IN_KPC; // puts it in units of km^2 / s^2

    let mut acceleration: [f64; 3] = [0.0; 3];
    let mut jerk: [f64; 3] = [0.0; 3];
    let mut snap: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        // put relative accel in km^2 / kpc s^2 so it has same units as other terms in snap sum
        relative_accel[i] *= KM_IN_KPC;

        // units of accel are naively km^2 / kpc s^2 but we want km / s^2
        acceleration[i] = GG * source.mass * displacement[i] / displacement_mag.powi(3);
        acceleration[i] /= KM_IN_KPC;
    
        // units of jerk are naively km^3 / kpc^2 s^2 but we want km / s^3
        jerk[i] = GG * &source.mass * ((relative_velocity[i] / displacement_mag.powi(3)) - 3.0 * rdotv * (displacement[i] / displacement_mag.powi(5)));
        jerk[i] /= KM_IN_KPC.powi(2);

        // units of snap are naively km^4 / kpc^3 s^3 but we want km / s^4
        snap[i] = GG * &source.mass *(
            (relative_accel[i] / displacement_mag.powi(3))
            - 6.0 * (rdotv * relative_velocity[i] / displacement_mag.powi(5))
            - 3.0 * (relative_velocity_squared + rdota) * (displacement[i] / displacement_mag.powi(5))
            + 15.0 * rdotv.powi(2) * (displacement[i] / displacement_mag.powi(7)) 
        );
        snap[i] /= KM_IN_KPC.powi(3);
    };

    (acceleration, jerk, snap)
}

pub fn calculate_energy(data: &Vec<Particle>) -> (f64, f64) {
    // units of energy are Msun km^2 / s^2
    let mut potential: f64 = 0.0;
    let mut kinetic: f64 = 0.0;
    for i in 0..data.len() {
        kinetic += kinetic_energy(&data[i]);

        for j in 0..data.len() {
            if i == j {
                continue;
            };
            potential += 0.5 * potential_of_binary_config(&data[i], &data[j]);
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

pub fn calculate_angular_momentum(data: &Vec<Particle>) -> [f64; 3] {
    // units of energy are Msun km kpc / s
    let mut angular_momentum: [f64; 3] = [0.0, 0.0, 0.0];
    for i in 0..data.len() {
        angular_momentum = add(&angular_momentum, &scalar_multiply(&data[i].mass, &cross_product(&data[i].position, &data[i].velocity)));
    };
    angular_momentum
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time_evol::compute_timestep;

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

        let dynamics = calculate_dynamics_due_to_one_body(&earth, &earth.acceleration.clone(), &sun);
        earth.acceleration = dynamics.0;
        earth.jerk = dynamics.1;
        earth.snap = dynamics.2;

        let timestep = compute_timestep(&earth, &0.1);
        println!("dt = {}", timestep);

        let expected_acceleration = 5.930e-6;
        let expected_jerk = 1.1805e-12;
        let expected_snap = 2.3484e-19;
        let expected_timestep = 3.55260e5;

        let accel_error = (magnitude(&dynamics.0)/expected_acceleration) - 1.0;
        let jerk_error = (magnitude(&dynamics.1)/expected_jerk) - 1.0;
        let snap_error = (magnitude(&dynamics.2)/expected_snap) - 1.0;
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
}