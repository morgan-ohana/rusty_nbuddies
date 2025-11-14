use crate::black_hole::BlackHole;

pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const KM_IN_KPC: f64 = 30856776000000000.0; // number of km in kpc

pub fn recalculate_dynamics_due_to_gravity(data: &mut Vec<BlackHole>) {
    for i in 0..data.len() {
        data[i].acceleration = [0.0; 3];
        for j in 0..data.len() {
            if i == j {
                continue;
            };
            let dynamical_variables = &calculate_dynamics_due_to_one_body(&data[i], &data[j]);
            data[i].acceleration = add(&data[i].acceleration, &dynamical_variables.0);
            data[i].jerk = add(&data[i].jerk, &dynamical_variables.1);
            data[i].snap = add(&data[i].snap, &dynamical_variables.2);
        };
    };
}

fn calculate_dynamics_due_to_one_body(target: &BlackHole, source: &BlackHole) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let displacement: [f64; 3] = subtract(&source.position, &target.position);
    let displacement_mag: f64 = magnitude(&displacement);
    
    let relative_velocity: [f64; 3] = subtract(&source.velocity, &target.velocity);
    let velocity_squared: f64 = dot_product(&relative_velocity, &relative_velocity);
    let rdotv: f64 = dot_product(&displacement, &relative_velocity);

    let mut relative_accel: [f64; 3] = subtract(&source.acceleration, &target.acceleration); // puts in units of km^2 / kpc s^2
    let rdota: f64 = dot_product(&displacement, &relative_accel) * KM_IN_KPC; // puts it in units of km^2 / s^2

    let mut acceleration: [f64; 3] = [0.0; 3];
    let mut jerk: [f64; 3] = [0.0; 3];
    let mut snap: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        // put relative accel in km^2 / kpc s^2 so it has same units as other terms in snap sum
        relative_accel[i] *= KM_IN_KPC;

        // units of accel are naively km^2 / kpc s^2 but we want km / s^2
        acceleration[i] = GG * source.mass * displacement[i] / (displacement_mag * displacement_mag * displacement_mag);
        acceleration[i] /= KM_IN_KPC;
    
        // units of jerk are naively km^3 / kpc^2 s^2 but we want km / s^3
        jerk[i] = GG * &source.mass * ((relative_velocity[i] / displacement_mag.powi(3)) - 3.0 * rdotv * (displacement[i] / displacement_mag.powi(5)));
        jerk[i] /= KM_IN_KPC.powi(2);

        // units of snap are naively km^4 / kpc^3 s^3 but we want km / s^4
        snap[i] = GG * &source.mass *(
            (relative_accel[i] / displacement_mag.powi(3))
            - 6.0 * rdotv * (relative_velocity[i] / displacement_mag.powi(5))
            - 3.0 * (velocity_squared + rdota) * (displacement[i] / displacement_mag.powi(5))
            + 15.0 * rdotv.powi(2) * (displacement[i] / displacement_mag.powi(7)) 
        );
        snap[i] /= KM_IN_KPC.powi(3);
    };

    (acceleration, jerk, snap)
}

pub fn calculate_energy(data: &Vec<BlackHole>) -> (f64, f64) {
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

fn kinetic_energy(bh: &BlackHole) -> f64 {
    0.5 * bh.mass * dot_product(&bh.velocity, &bh.velocity)
}

fn potential_of_binary_config(bh1: &BlackHole, bh2: &BlackHole) -> f64 {
    -1.0 * GG * bh1.mass * bh2.mass / magnitude(&subtract(&bh2.position, &bh1.position))
}

pub fn calculate_angular_momentum(data: &Vec<BlackHole>) -> [f64; 3] {
    // units of energy are Msun km kpc / s
    let mut angular_momentum: [f64; 3] = [0.0, 0.0, 0.0];
    for i in 0..data.len() {
        angular_momentum = add(&angular_momentum, &scalar_multiply(&data[i].mass, &cross_product(&data[i].position, &data[i].velocity)));
    };
    angular_momentum
}

// Vector Utilities

fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let mut c: [f64; 3] = [0.0, 0.0, 0.0];
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[2];
    c
}

fn magnitude(vec: &[f64; 3]) -> f64 {
    f64::sqrt(dot_product(vec, vec))
}

fn add(vec1: &[f64; 3], vec2: &[f64; 3]) -> [f64; 3] {
    let mut sum: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        sum[i] = vec1[i] + vec2[i];
    }
    sum
}

fn subtract(vec1: &[f64; 3], vec2: &[f64; 3]) -> [f64; 3] {
    let mut sum: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        sum[i] = vec1[i] - vec2[i];
    }
    sum
}

fn scalar_multiply(scalar: &f64, vec: &[f64; 3]) -> [f64; 3] {
    let mut output = vec.clone();
    for i in 0..3 {
        output[i] = scalar * vec[i];
    }
    output
}
