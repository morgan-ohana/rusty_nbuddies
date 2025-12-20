use std::f64::consts::PI;
use std::f64::consts::SQRT_2;
use crate::forces::GG;
use crate::init_conds::SPACIAL_GRID_NUM;
use crate::init_conds::VELOCITY_GRID_NUM;

use crate::plotting::plot_check_function;
use crate::plotting::plot_function;

const ENERGY_GRID_NUM: usize = 10000;
const OFFSET: f64 = 1e-10; //for when something really just shouldn't be zero

pub fn compute_phase_space_density(rho_points: &Vec<f64>, v_points: &Vec<f64>, potential_points: &Vec<f64>, cuspy: &bool) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    plot_check_function(&potential_points, &|v: f64| -> anyhow::Result<f64> {Ok(-3.0 * v.powi(5) / (4.0 * PI * GG.powi(5) * 1e8_f64.powi(4)))}, &rho_points, &"rho_vs_V_test.png", &"drho/dV", &"V", &"rho")?;

    let d2rho_dpotential2 = differentiate(&differentiate(&rho_points, &potential_points, cuspy, true), &potential_points, cuspy, false);
    
    let energy_min: f64 = potential_points[0]*1.001; // .1% buffer to prevent any points having exactly the minimum energy;
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

    let mut f_of_energy: Vec<f64> = vec![0.0; ENERGY_GRID_NUM];

    //values at max r and v left at 0 (presumably those regions are unoccupied anyway)
    let mut integral_range_exhausted = false;
    for j in 0..ENERGY_GRID_NUM {
        for k in 2..(d2rho_dpotential2.len()) {
            //Index backwards from E_max at r_max down to E at r, start at len - 2 since int is computed by trapezoids between i and i+1
            let i = d2rho_dpotential2.len() - k;
            //Note d2r_d2V is one point shorter on both ends so V[i+1] is at the same place as d2rdV2[i]
            let l = i + 1;

            //enforce integral lower bound
            if potential_points[i+1] <= energy_points[j] {
                if k == 1 && !integral_range_exhausted {
                    println!("integral range exhausted at {}", j);
                    integral_range_exhausted = true;
                }
                break
            }

            let delta_potential = potential_points[l+1] - potential_points[l];
            let integrand_term_1 = d2rho_dpotential2[i] / (potential_points[l] - energy_points[j]).sqrt();
            let integrand_term_2 = d2rho_dpotential2[i+1] / (potential_points[l+1] - energy_points[j]).sqrt();
            
            //trapezoidal sum
            f_of_energy[j] += 0.5 * (integrand_term_1 + integrand_term_2) * delta_potential;
        }
        //prefactor
        f_of_energy[j] /= PI.powi(2) * 2.0 * SQRT_2;
        
        if f_of_energy[j] < 0.0 {
            panic!("Negative distribution function, Eddington Inversion has failed at j={}", j);
        }
    }
        

    plot_check_function(&energy_points, &|e: f64| -> anyhow::Result<f64> {Ok((24.0*2.0_f64.sqrt()/(7.0*PI.powi(3))) * (GG.powi(-5))*(1e8_f64).powi(-4) * (-e).powi(7).sqrt())}, &f_of_energy, &"Eddington_Inversion_Phasespace_Distribution_test.png", &"Phase Space Distribution From Eddington Inverter", &"Energy (M_sun km^2/s^2)", &"f (M_sun s^3/kpc^3 km^3)")?;

    let mut f: Vec<Vec<f64>> = vec![vec![0.0; VELOCITY_GRID_NUM]; SPACIAL_GRID_NUM];

    for i in 0..SPACIAL_GRID_NUM {
        for j in 0..VELOCITY_GRID_NUM {
            let energy: f64 = 0.5 * v_points[j] * v_points[j] + potential_points[i];
            
            //ignore unbound regions, f=0 there anyway
            if energy > 0.0 {
                continue
            }
            
            //binary search for neighborhood
            let mut low = 0;
            let mut high = ENERGY_GRID_NUM - 1;
            
            while high - low > 1 {
                let mid = (low + high) / 2;

                //check if above or below
                match energy > energy_points[mid] {
                    true => low = mid,
                    false => high = mid
                };
            }
            //low and high now bracket the energy value
            
            //linear interpolation over interval
            let t = (energy - energy_points[low]) / (energy_points[high] - energy_points[low]);
            f[i][j] = f_of_energy[low] * (1.0 - t) + f_of_energy[high] * t;
            
        }
    }

    Ok(f)
}

pub fn generate_rho_points<T: Fn(f64) -> anyhow::Result<f64>>(rho: &T, r_points: &Vec<f64>) -> anyhow::Result<Vec<f64>> {
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

fn differentiate(y: &Vec<f64>, x: &Vec<f64>, cuspy: &bool, forward: bool) -> Vec<f64> {
    //Assumes y and x start at same point, they do not need to end at same point. That is reverse(forward) is safe, forward(reverse) is not
    let mut dy_dx: Vec<f64> = vec![0.0; y.len() - 1];
    let start = match forward {
        true => 0,
        false => 1
    };

    let end = match forward {
        true => y.len() - 1,
        false => y.len()
    };

    for i in start..end {
        let j = match forward {
            true => i,
            false => i-1
        };

        dy_dx[j] = match cuspy {
            true => {
                // using dy/dx = (y + offset) dlog(y + OFFSET)/dx with offset to avoid logs of 0 anywhere
                (0.5*(y[j+1] + y[j])) * ((y[j+1] + OFFSET).ln() - (y[j] + OFFSET).ln()) / (x[j+1] - x[j])
            }
            false => {
                (y[j+1] - y[j]) / (x[j+1] - x[j])
            }
        }
    }
    dy_dx
}