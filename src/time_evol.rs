use crate::black_hole::BlackHole;
use crate::forces::KM_PER_KPC;

pub fn update_positions(data: &mut Vec<BlackHole>, previous_data: &Vec<BlackHole>, delta_t: &f64) {
    for n in 0..data.len() {
        for i in 0..3 {
            // velocity/accel in km/s system but want pos in kpc
            data[n].position[i] += (previous_data[n].velocity[i] * delta_t + 0.5 * previous_data[n].acceleration[i] * delta_t * delta_t)/KM_PER_KPC;
        }
    }
}

pub fn update_velocities(data: &mut Vec<BlackHole>, previous_data: &Vec<BlackHole>, delta_t: &f64) {
    for n in 0..data.len() {
        for i in 0..3 {
            data[n].velocity[i] += 0.5 * (data[n].acceleration[i] + previous_data[n].acceleration[i]) * delta_t;
        }
    }
}
