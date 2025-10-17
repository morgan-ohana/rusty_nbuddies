use num_traits::pow;
use std::f64::consts::PI;
use crate::forces::GG;

const SPACIAL_GRID_NUM: usize = 1000;
const VELOCITY_GRID_NUM: usize = 1000;
const ENERGY_GRID_NUM: usize = 1000;
const LOG_OFFSET: f64 = 1e-10;

fn compute_phase_space_density(rho: fn(f64) -> f64, r_cutoff: &f64, v_cutoff: &f64, mass: &f64, cuspy: bool) -> anyhow::Result<[[f64; VELOCITY_GRID_NUM]; SPACIAL_GRID_NUM]> {
    let mut r_points: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];
    for i in 0..SPACIAL_GRID_NUM {
        r_points[i] = (i as f64) * r_cutoff / ((SPACIAL_GRID_NUM - 1) as f64);
    }

    let rho_points = generate_rho_points(rho, &r_points);
    let potential_points = generate_potential_points(&rho_points, &r_points);
    let drho_dV = generate_drho_dV(&rho_points, &potential_points, cuspy);

    let energy_cutoff: f64 = 0.5 * mass * v_cutoff * v_cutoff + potential_points[SPACIAL_GRID_NUM - 1]; //max potential occurs at largest radius

    let mut energy_points: [f64; ENERGY_GRID_NUM] = [0.0; ENERGY_GRID_NUM];
    for j in 0..ENERGY_GRID_NUM {
        energy_points[j] = potential_points[0] + (j as f64) * (energy_cutoff - potential_points[0]) / ((ENERGY_GRID_NUM - 1) as f64); //min energy is min potential at r=0
    }
    let mut integral: [f64; ENERGY_GRID_NUM] = [0.0; ENERGY_GRID_NUM];

    //values at max r and v left at 0 (presumably those regions are unoccupied anyway)
    for j in 0..ENERGY_GRID_NUM {
        for k in 1..SPACIAL_GRID_NUM {
            //Index backwards from E_max at r_max down to E at r
            let i = SPACIAL_GRID_NUM - 1 - k;
        
            let delta_V = potential_points[i+1] - potential_points[i];
            let integrand = drho_dV[i] / (potential_points[i] - energy_points[j]).sqrt();
            
            integral[j] += integrand * delta_V;
        }
    }

    let mut f_of_energy: [f64; ENERGY_GRID_NUM] = [0.0; ENERGY_GRID_NUM];

    for j in 0..(ENERGY_GRID_NUM - 1) {
        f_of_energy[j] = (integral[j+1] - integral[j]) / (energy_points[j+1] - energy_points[j]);
        f_of_energy[j] /= 8.0_f64.sqrt() * PI * PI;
        if f_of_energy[j] < 0.0 {
            panic!("Negative distribution function, Eddington Inversion has failed");
        }
    }
    
    let mut v_points: [f64; VELOCITY_GRID_NUM] = [0.0; VELOCITY_GRID_NUM];
    for j in 0..VELOCITY_GRID_NUM {
        v_points[j] = (j as f64) * v_cutoff / ((VELOCITY_GRID_NUM - 1) as f64);
    }

    let mut f: [[f64; ENERGY_GRID_NUM]; SPACIAL_GRID_NUM] = [[0.0; ENERGY_GRID_NUM]; SPACIAL_GRID_NUM];

    for i in 0..SPACIAL_GRID_NUM {
        for j in 0..VELOCITY_GRID_NUM {
            let energy: f64 = 0.5 * mass * v_points[j] * v_points[j] + potential_points[i];
            let mut e_arg: usize = 0;
            for k in 0..ENERGY_GRID_NUM {
                if (energy_points[k] - energy).abs() < (energy_points[e_arg] - energy).abs() {
                    e_arg = k;
                }
            }
            f[i][j] = f_of_energy[e_arg];
        }
    }

    Ok(f)
}

fn generate_rho_points(rho: fn(f64) -> f64, r_points: &[f64; SPACIAL_GRID_NUM]) -> [f64; SPACIAL_GRID_NUM] {
    let mut rho_points: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];
    for i in 0..SPACIAL_GRID_NUM {
        rho_points[i] = rho(r_points[i]);
    }

    rho_points
}

fn generate_potential_points(rho_points: &[f64; SPACIAL_GRID_NUM], r_points: &[f64; SPACIAL_GRID_NUM]) -> [f64; SPACIAL_GRID_NUM] {
    let mut mass_enclosed: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];
    
    mass_enclosed[0] = rho_points[0] * 4.0 * PI * pow(0.5*r_points[1], 3)/3.0;
    for i in 1..(SPACIAL_GRID_NUM - 1) {
        mass_enclosed[i] = mass_enclosed[i-1] + rho_points[i] * 4.0 * PI * (pow(0.5*(r_points[i+1] + r_points[i]), 3) - pow(0.5*(r_points[i] + r_points[i-1]), 3)) / 3.0;
    }
    // Assumes linear spacing or r for last step, shouldn't matter since rho should basically be zero out here anyway
    mass_enclosed[SPACIAL_GRID_NUM - 1] = mass_enclosed[SPACIAL_GRID_NUM - 2] + rho_points[SPACIAL_GRID_NUM - 1] * 4.0 * PI * (pow(1.5*r_points[SPACIAL_GRID_NUM - 1] - 0.5*r_points[SPACIAL_GRID_NUM - 2], 3) - pow(0.5*(r_points[SPACIAL_GRID_NUM - 1] + r_points[SPACIAL_GRID_NUM-2]), 3)) / 3.0;

    let mut potential_points: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];

    potential_points[SPACIAL_GRID_NUM - 1] = -1.0 * (r_points[SPACIAL_GRID_NUM - 1] - r_points[SPACIAL_GRID_NUM - 2]) * GG * mass_enclosed[SPACIAL_GRID_NUM - 1] / (r_points[SPACIAL_GRID_NUM-1] * r_points[SPACIAL_GRID_NUM-1]);
    for i in 1..SPACIAL_GRID_NUM {
        let j = SPACIAL_GRID_NUM - 1 - i;
        let delta_r = r_points[j+1] - r_points[j-1];
        let force = - 1.0 * GG * mass_enclosed[j] / (r_points[j]*r_points[j]);
        potential_points[j] = potential_points[j+1] + delta_r * force;
    }

    potential_points
}

fn generate_drho_dV(rho_points: &[f64; SPACIAL_GRID_NUM], potential_points: &[f64; SPACIAL_GRID_NUM], cuspy: bool) -> [f64; SPACIAL_GRID_NUM] {
    let mut drho_dV: [f64; SPACIAL_GRID_NUM] = [0.0; SPACIAL_GRID_NUM];

    for i in 0..(SPACIAL_GRID_NUM - 1) {
        drho_dV[i] = match cuspy {
            true => {
                // using drho/dV = log(rho + LOG_OFFSET) dlog(rho + LOG_OFFSET)/dV with offset to avoid logs of 0 anywhere
                (0.5*(rho_points[i+1] + rho_points[i]) + LOG_OFFSET).ln() * ((rho_points[i+1] + LOG_OFFSET).ln() - (rho_points[i] + LOG_OFFSET).ln()) / (potential_points[i+1] - potential_points[i])
            }
            false => {
                (rho_points[i+1] - rho_points[i]) / (potential_points[i+1] - potential_points[i])
            }
        }
    }
    //final value left at 0, these large r values shouldn't matter anyway.

    drho_dV
}