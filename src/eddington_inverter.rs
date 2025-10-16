use num_traits::pow;
use std::f64::consts::PI;
use crate::forces::GG;

const SPACIAL_GRID_NUM: usize = 1000;

fn generate_rho_points(rho: fn(f64) -> f64, r_points: &[f64; SPACIAL_GRID_NUM]) -> [f64; SPACIAL_GRID_NUM] {
    let mut rho_points: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];
    for i in 0..SPACIAL_GRID_NUM {
        rho_points[i] = rho(r_points[i]);
    }

    rho_points
}

fn generate_V_points(rho_points: &[f64; SPACIAL_GRID_NUM], r_points: &[f64; SPACIAL_GRID_NUM]) -> [f64; SPACIAL_GRID_NUM] {
    let mut M_enclosed: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];
    
    M_enclosed[0] = rho_points[0] * 4.0 * PI * pow(0.5*r_points[1], 3)/3.0;
    for i in 1..(SPACIAL_GRID_NUM - 1) {
        M_enclosed[i] = M_enclosed[i-1] + rho_points[i] * 4.0 * PI * (pow(0.5*(r_points[i+1] + r_points[i]), 3) - pow(0.5*(r_points[i] + r_points[i-1]), 3)) / 3.0;
    }
    // Assumes linear spacing or r for last step, shouldn't matter since rho should basically be zero out here anyway
    M_enclosed[SPACIAL_GRID_NUM - 1] = M_enclosed[SPACIAL_GRID_NUM - 2] + rho_points[SPACIAL_GRID_NUM - 1] * 4.0 * PI * (pow(1.5*r_points[SPACIAL_GRID_NUM - 1] - 0.5*r_points[SPACIAL_GRID_NUM - 2], 3) - pow(0.5*(r_points[SPACIAL_GRID_NUM - 1] + r_points[SPACIAL_GRID_NUM-2]), 3)) / 3.0;

    let mut V_points: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];

    V_points[SPACIAL_GRID_NUM - 1] = -1.0 * (r_points[SPACIAL_GRID_NUM - 1] - r_points[SPACIAL_GRID_NUM - 2]) * GG * M_enclosed[SPACIAL_GRID_NUM - 1] / (r_points[SPACIAL_GRID_NUM-1] * r_points[SPACIAL_GRID_NUM-1]);
    for i in 1..SPACIAL_GRID_NUM {
        let j = SPACIAL_GRID_NUM - 1 - i;
        let delta_r = 0.5*(r_points[j+1] - r_points[j-1]);
        let force = - 1.0 * GG * M_enclosed[j] / (r_points[j]*r_points[j]);
        V_points[j] = V_points[j+1] + delta_r * force;
    }

    V_points
}

fn generate_drho_dV(rho_points: &[f64; SPACIAL_GRID_NUM], V_points: &[f64; SPACIAL_GRID_NUM]) -> [f64; SPACIAL_GRID_NUM] {
    let mut drho_dV: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];

    for i in 0..(SPACIAL_GRID_NUM - 1) {
        drho_dV[i] = (rho_points[i+1] - rho_points[i]) / (V_points[i+1] - V_points[i])
    }
    //final value left at 0, these large r values shouldn't matter anyway.

    drho_dV
}