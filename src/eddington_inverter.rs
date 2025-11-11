use num_traits::pow;
use std::f64::consts::PI;
use crate::forces::GG;
use crate::init_conds::SPACIAL_GRID_NUM;
use crate::init_conds::VELOCITY_GRID_NUM;

use crate::plotting::plot_func;

const ENERGY_GRID_NUM: usize = 10000;
const OFFSET: f64 = 1e-10; //for when something really just shouldn't be zero

pub fn compute_phase_space_density(rho_points: &Vec<f64>, v_points: &Vec<f64>, potential_points: &Vec<f64>, cuspy: &bool) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    
    let drho_dpotential = generate_drho_dpotential(&rho_points, &potential_points, cuspy);
    
    let energy_min: f64 = potential_points[0];
    let energy_max: f64 = match cuspy {
        true => energy_min*1e-3,
        false => 0.0
    };

    let log_space_energy = false; //cuspy;
    let mut energy_points: Vec<f64> = vec![0.0; ENERGY_GRID_NUM];
    for j in 0..ENERGY_GRID_NUM {
        energy_points[j] = match log_space_energy {
            true => {
                -1.0 * ((-1.0*energy_min).ln() + (j as f64) * ((-1.0 * energy_max).ln() - (-1.0 * energy_min).ln()) / ((ENERGY_GRID_NUM - 1) as f64)).exp()
            }
            false => {
                energy_min + (j as f64) * (0.0 - energy_min) / ((ENERGY_GRID_NUM - 1) as f64)
            }
        };
    }

    let mut integral: Vec<f64> = vec![0.0; ENERGY_GRID_NUM];

    //values at max r and v left at 0 (presumably those regions are unoccupied anyway)
    let mut integral_range_exhausted = false;
    for j in 0..ENERGY_GRID_NUM {
        for k in 1..SPACIAL_GRID_NUM {
            //Index backwards from E_max at r_max down to E at r
            let i = SPACIAL_GRID_NUM - 1 - k;

            if potential_points[i] <= energy_points[j] {
                //enforce integral lower bound
                if k == 1 && !integral_range_exhausted {
                    println!("integral range exhausted at {}", j);
                    integral_range_exhausted = true;
                }
                break
            }

            let delta_potential = potential_points[i+1] - potential_points[i];
            let integrand = drho_dpotential[i] / (potential_points[i] - energy_points[j]).sqrt();
            
            integral[j] += integrand * delta_potential;
        }
    }
    
    denoise(&mut integral);
    
    let mut f_of_energy: Vec<f64> = vec![0.0; ENERGY_GRID_NUM];

    for j in 0..(ENERGY_GRID_NUM - 1) {
        f_of_energy[j] = match cuspy {
            true => (0.5 * (integral[j+1] + integral[j])) * ((-1.0*integral[j+1] + OFFSET).ln() - (-1.0*integral[j] + OFFSET).ln()) / (energy_points[j+1] - energy_points[j]),
            false => (integral[j+1] - integral[j]) / (energy_points[j+1] - energy_points[j])
        };
        f_of_energy[j] /= 8.0_f64.sqrt() * PI * PI;
        if f_of_energy[j] < 0.0 {
            panic!("Negative distribution function, Eddington Inversion has failed at j={}", j);
        }
    }
    
    let mut f: Vec<Vec<f64>> = vec![vec![0.0; VELOCITY_GRID_NUM]; SPACIAL_GRID_NUM];

    for i in 0..SPACIAL_GRID_NUM {
        for j in 0..VELOCITY_GRID_NUM {
            let energy: f64 = 0.5 * v_points[j] * v_points[j] + potential_points[i];
            
            //ignore unbound regions, f=0 there anyway
            if energy > 0.0 {
                continue
            }
            
            //binary search for neighborhood
            let mut search_space = ENERGY_GRID_NUM as f64;
            let mut e_arg: usize = (search_space / 2.0).round() as usize;
            
            while search_space >= 1.0 && e_arg != 0 && e_arg != ENERGY_GRID_NUM - 1 {
                //shrink search
                search_space /= 2.0;

                //enforce bounds on idx
                if search_space > e_arg as f64 {
                    search_space = e_arg as f64
                }
                if search_space > (ENERGY_GRID_NUM - 1 - e_arg) as f64 {
                    search_space = (ENERGY_GRID_NUM - 1 - e_arg) as f64
                }

                //check if above or below
                match energy > energy_points[e_arg] {
                    true => e_arg += (search_space / 2.0).round() as usize,
                    false => e_arg -= (search_space /2.0).round() as usize
                };
            }
            
            //check if at end and skip interpolation:
            if e_arg == ENERGY_GRID_NUM - 1 {
                f[i][j] = energy_points[e_arg];
                continue
            }

            //identity interval where E lies
            let mut found_interval = false;
            const NEIGHBORHOOD_WIDTH: usize= 1;
            for k in (e_arg - NEIGHBORHOOD_WIDTH)..(e_arg + NEIGHBORHOOD_WIDTH + 1) {
                if (energy_points[k+1] - energy) > 0.0 && (energy - energy_points[k]) >= 0.0 {
                    e_arg = k;
                    found_interval = true;
                    break
                }
            }
            if !found_interval {
                panic!("Failed to find correct interval for E at i = {}, j = {} with energy of {} in the neighborhood of {} with e_arg = {}", i, j, energy, energy_points[e_arg], e_arg)
            }

            //linear interpolation over interval
            let t = (energy - energy_points[e_arg]) / (energy_points[e_arg + 1] - energy_points[e_arg]);
            f[i][j] = f_of_energy[e_arg] * (1.0 - t) + f_of_energy[e_arg + 1] * t;
        }
    }

    Ok(f)
}

fn denoise(points: &mut Vec<f64>) {
    //Removes negative spikes patching them over with linear segments
    let mut j;
    for i in 0..points.len() {
        j = 0;
        while j < (points.len() - i - 1) && points[i+j+1] < points[i] {
            //println!("i={}, j={}", i, j);
            j += 1;
        }
        for k in 1..(j+1) {
            points[i+k] = points[i] + (points[i+j+1] - points[i])*((k as f64)/((j+1) as f64));
        }
    }
}

pub fn generate_rho_points<T: Fn(f64) -> anyhow::Result<f64>>(rho: T, r_points: &Vec<f64>) -> anyhow::Result<Vec<f64>> {
    let mut rho_points: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    for i in 0..SPACIAL_GRID_NUM {
        rho_points[i] = rho(r_points[i])?;
    }

    Ok(rho_points)
}

pub fn generate_potential_points(rho_points: &Vec<f64>, r_points: &Vec<f64>) -> Vec<f64> {
    //Uses direct formula V = -4pi G (1/r int_0^r rho r'^2 dr' + int_r^inf rho r' dr')
    let mut potential_points: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];

    for i in 0..SPACIAL_GRID_NUM {
        //V = 0
        
        //linear approx for center
        let rho_at_center: f64 = rho_points[0] - r_points[0] * ((rho_points[1] - rho_points[0])/(r_points[1] - r_points[0]));
        potential_points[i] += 0.5 * (rho_points[0] + rho_at_center) * r_points[0].powi(3) / 3.0;
        for j in 1..(i+1) {
            potential_points[i] += 0.5 * (rho_points[j] + rho_points[j-1]) * (r_points[j].powi(3) - r_points[j-1].powi(3)) / 3.0;
        }
        //V = int_0^r rho r'^2 dr'
        
        potential_points[i] /= r_points[i] + OFFSET;
        //V = (1/r) int_0^r rho r'^2 dr'
        
        for j in i..(SPACIAL_GRID_NUM - 1) {
            potential_points[i] += 0.5 * (rho_points[j+1] + rho_points[j]) * (r_points[j+1].powi(2) - r_points[j].powi(2)) / 2.0;
        }
        //V = (1/r) int_0^r rho r'^2 dr' + int_r^inf rho r' dr'

        potential_points[i] *= -4.0 * PI * GG;
        //V = -4pi G ((1/r) int_0^r rho r'^2 dr' + int_r^inf rho r' dr')
    }

    potential_points
}

fn generate_drho_dpotential(rho_points: &Vec<f64>, potential_points: &Vec<f64>, cuspy: &bool) -> Vec<f64> {
    let mut drho_dpotential: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];

    for i in 0..(SPACIAL_GRID_NUM - 1) {
        drho_dpotential[i] = match cuspy {
            true => {
                // using drho/dV = log(rho + OFFSET) dlog(rho + OFFSET)/dV with offset to avoid logs of 0 anywhere
                (0.5*(rho_points[i+1] + rho_points[i])) * ((rho_points[i+1] + OFFSET).ln() - (rho_points[i] + OFFSET).ln()) / (potential_points[i+1] - potential_points[i])
            }
            false => {
                (rho_points[i+1] - rho_points[i]) / (potential_points[i+1] - potential_points[i])
            }
        }
    }
    //final value left at 0, these large r values shouldn't matter anyway.

    drho_dpotential
}