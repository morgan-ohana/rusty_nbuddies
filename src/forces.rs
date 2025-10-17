use crate::black_hole::BlackHole;

pub const GG: f64 = 4.301e-6; // Newton constant km^2 kpc / Msun s^2
pub const KM_PER_KPC: f64 = 30856776000000000.0; // number of km in kpc

pub fn recalculate_acceleration_due_to_gravity(data: &mut Vec<BlackHole>) {
    for i in 0..data.len() {
        data[i].acceleration = [0.0; 3];
        for j in 0..data.len() {
            if i == j {
                continue;
            };
            data[i].acceleration = add(&data[i].acceleration, &calculate_acceleration_due_to_one_body(&data[i], &data[j]));
        };
    };
}

fn calculate_acceleration_due_to_one_body(target: &BlackHole, source: &BlackHole) -> [f64; 3] {
    let displacement: [f64; 3] = subtract(&source.position, &target.position);
    let displacement_mag: f64 = magnitude(&displacement);
    let mut acceleration: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        // units of accel are naively km^2 / kpc s^2 but we want km / s^2
        acceleration[i] = GG * source.mass * displacement[i] / (displacement_mag * displacement_mag * displacement_mag);
        acceleration[i] /= KM_PER_KPC;
    };
    acceleration
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
