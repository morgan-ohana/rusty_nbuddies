use crate::particle::Particle;
use crate::forces::KM_IN_KPC;
use crate::vectors::magnitude;

pub fn update_positions(data: &mut Vec<Particle>, delta_t: &f64) {
    for n in 0..data.len() {
        for i in 0..3 {
            // velocity/accel in km/s system but want pos in kpc
            data[n].position[i] += (data[n].velocity[i] * delta_t) /KM_IN_KPC;
        }
    }
}

pub fn update_velocities(data: &mut Vec<Particle>, delta_t: &f64) {
    for n in 0..data.len() {
        for i in 0..3 {
            data[n].velocity[i] += data[n].acceleration[i] * delta_t;
        }
    }
}

pub fn compute_timestep(particle: &Particle, eta: &f64) -> f64 {
    let accel = magnitude(&particle.acceleration);
    let jerk = magnitude(&particle.jerk);
    let snap = magnitude(&particle.snap);
    // println!("a = {}", accel);
    // println!("j = {}", jerk);
    // println!("s = {}", snap);

    eta / ((jerk/accel).powi(2) + (snap/accel)).sqrt()
}