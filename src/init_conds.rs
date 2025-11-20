use std::f64::consts::PI;
use plotters::prelude::*;
use num_traits::pow;

use crate::black_hole::BlackHole;
use crate::forces::GG;
use crate::forces::KM_IN_KPC;
use crate::eddington_inverter::*;
use crate::plotting::plot_func;

pub const SPACIAL_GRID_NUM: usize = 10000;
pub const VELOCITY_GRID_NUM: usize = 10000;

pub fn abg_profile_init_conds(alpha: &f64, beta: &f64, gamma: &f64, r_s: &f64, rho_s: &f64, r_cutoff: &f64, particle_num: &usize) -> Result<(), Box<dyn std::error::Error>> {

    let r_max = 10.0*r_cutoff; // arbitrary factor

    let v_max = 400.0;

    let cuspy: bool = gamma > &0.0;

    let abg_profile = |r: f64| -> anyhow::Result<f64> {
        if *gamma > 3.0 {
            panic!("Total mass diverges at small r");
        }
        if *beta > 3.0 || r < *r_cutoff {
            Ok(rho_s / ((r / r_s).powf(*gamma) * (1.0 + (r / r_s).powf(*alpha)).powf((*beta - *gamma) / *alpha)))
        } else {
            let rho_cutoff = rho_s / ((r_cutoff / r_s).powf(*gamma) * (1.0 + (r_cutoff / r_s).powf(*alpha)).powf((*beta - *gamma) / *alpha));
            let r_decay = 0.3 * r_cutoff;
            let delta = (r_cutoff / r_decay) - ((*gamma + *beta * (r_cutoff / r_s).powf(*alpha)) / (1.0 + (r_cutoff / r_s).powf(*alpha)));
            Ok(rho_cutoff * (r / r_cutoff).powf(delta) * ((r_cutoff - r) / r_decay).exp())
        }
    };

    generate_init_conds_from_rho(abg_profile, &r_max, &particle_num, &cuspy)?;

    Ok(())
}

pub fn plummer_init_conds(r_s: &f64, total_mass: &f64, particle_num: &usize) -> Result<(), Box<dyn std::error::Error>> {
    let r_max = 3.0*r_s; //3.0*r_s;
    let v_max = 400.0;

    let cuspy = false;

    let plummer_profile = |r: f64| -> anyhow::Result<f64> {
        Ok(3.0*(*total_mass)*(*r_s)*(*r_s)/(4.0*PI*((*r_s)*(*r_s) + r*r).powf(2.5)))
    };

    generate_init_conds_from_rho(plummer_profile, &r_max, &particle_num, &cuspy)?;

    Ok(())
}

fn generate_init_conds_from_rho<T: Fn(f64) -> anyhow::Result<f64>>(rho: T, r_max: &f64, particle_num: &usize, cuspy: &bool) -> Result<(), Box<dyn std::error::Error>> {
    let mut r_points: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    for i in 0..SPACIAL_GRID_NUM {
        r_points[i] = match cuspy {
            true => {
                let r_min = r_max*1e-4;
                (r_min.ln() + (i as f64) * (r_max.ln() - r_min.ln()) / ((SPACIAL_GRID_NUM - 1) as f64)).exp()
            }
            false => {
                (i as f64) * r_max / ((SPACIAL_GRID_NUM - 1) as f64)
            }
        }
    }

    let rho_points = generate_rho_points(rho, &r_points)?;

    let potential_points = generate_potential_points(&rho_points, &r_points);
    let v_max = (-2.0 * potential_points[0]).sqrt();

    let mut v_points: Vec<f64> = vec![0.0; VELOCITY_GRID_NUM];
    for j in 0..VELOCITY_GRID_NUM {
        v_points[j] = (j as f64) * v_max / ((VELOCITY_GRID_NUM - 1) as f64);
    }
    println!("v_max={}", v_points[VELOCITY_GRID_NUM - 1]);

    let mut p_of_r: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    let mut total_mass: f64 = 0.0;

    p_of_r[0] = rho_points[0] * 4.0 * PI * pow(0.5*r_points[1], 3)/3.0;
    total_mass += p_of_r[0];
    for i in 1..(SPACIAL_GRID_NUM - 1) {
        p_of_r[i] = rho_points[i] * 4.0 * PI * (pow(0.5*(r_points[i+1] + r_points[i]), 3) - pow(0.5*(r_points[i] + r_points[i-1]), 3)) / 3.0;
        total_mass += p_of_r[i];
    }
    // Assumes linear spacing of r for last step, shouldn't matter since rho should basically be zero out here anyway
    p_of_r[SPACIAL_GRID_NUM - 1] = rho_points[SPACIAL_GRID_NUM - 1] * 4.0 * PI * (pow(1.5*r_points[SPACIAL_GRID_NUM - 1] - 0.5*r_points[SPACIAL_GRID_NUM - 2], 3) - pow(0.5*(r_points[SPACIAL_GRID_NUM - 1] + r_points[SPACIAL_GRID_NUM-2]), 3)) / 3.0;
    total_mass += p_of_r[SPACIAL_GRID_NUM - 1];

    //Normalizing
    for i in 0..SPACIAL_GRID_NUM {
        p_of_r[i] /= total_mass;
    }

    let f: Vec<Vec<f64>> = compute_phase_space_density(&rho_points, &v_points, &potential_points, &cuspy)?;

    let mut p_of_v_given_r: Vec<Vec<f64>> = vec![vec![0.0; VELOCITY_GRID_NUM]; SPACIAL_GRID_NUM];
    let mut normalization_check: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    
    let mut norm_drop: bool = false;
    for i in 0..SPACIAL_GRID_NUM {
        let v_escape = (-2.0 * potential_points[i]).sqrt();
        //println!("v_esc = {}, v_max = {}", v_escape, v_max);
        let v_esc_arg = ((v_escape/v_max) * VELOCITY_GRID_NUM as f64) as usize - 1;
        let mut non_phys_prob = 0.0;

        p_of_v_given_r[i][0] = (4.0 * PI / 3.0) * (f[i][0] / rho_points[i]) * (0.5*v_points[1]).powi(3);
        normalization_check[i] += p_of_v_given_r[i][0];
        
        for j in 1..v_esc_arg {
            
            p_of_v_given_r[i][j] = (4.0 * PI / 3.0) * (f[i][j] / rho_points[i]) * ((0.5*(v_points[j+1] + v_points[j])).powi(3) - (0.5*(v_points[j] + v_points[j-1])).powi(3));
            normalization_check[i] += p_of_v_given_r[i][j];
            
            if j >= v_esc_arg {
                panic!("bounds not working");
            }
        }

        println!("Normalization check: {}", normalization_check[i]);
        //println!("p(v>v_esc|r) = {}", non_phys_prob);
        if normalization_check[i] < 1.0 && !norm_drop {
            norm_drop = true;
            println!("Possible V_max bounds issue at r = {} with p(v|r)={} at the end", r_points[i], p_of_v_given_r[i][VELOCITY_GRID_NUM - 1])
        }

        //renormalizing
        for j in 0..VELOCITY_GRID_NUM {
            //p_of_v_given_r[i][j] /= normalization_check[i];
        }
    }

    plot_func(&r_points, &compute_velocity_dispersion(&v_points, &p_of_v_given_r), &"test_v_disp.png")?;
    plot_rho_recovery_test(&r_points, &rho_points, &recover_rho(&f, &v_points), &"Eddington_Inversion_convergence_test.png")?;

    Ok(())
}

pub fn binary_init_conds(seperation: f64, angular_momentum: f64, m1: f64, m2: f64) -> Vec<BlackHole> {
    // 0 = p1 + p2
    // l = p1 r1 + p2 r2 = p1 (r1 - r2) = p1 sep 
    // v1 = l/m1 sep
    // v2 = -l/m2 sep
    let bh_1 = BlackHole {
        mass: m1,
        position: [0.5*seperation, 0.0, 0.0],
        velocity: [0.0, angular_momentum / (m1 * seperation), 0.0],
        acceleration: [0.0; 3],
        jerk: [0.0; 3],
        snap: [0.0; 3],
    };
    let bh_2 = BlackHole {
        mass: m2,
        position: [-0.5*seperation, 0.0, 0.0],
        velocity: [0.0, - angular_momentum / (m2 * seperation), 0.0],
        acceleration: [0.0; 3],
        jerk: [0.0; 3],
        snap: [0.0; 3],
    };
    vec![bh_1, bh_2]
}

pub fn binary_circular_init_conds(seperation: f64, m1: f64, m2: f64) -> Vec<BlackHole> {
    // 0 = m1*x + m2 (x - sep) = - m2 sep + (m1 + m2)x => x = m2 sep/ (m1 + m2)
    let offset = m2 * seperation / (m1+ m2);
    // G (m1m2/mu) / sep^2 = G (m1 + m2) / sep ^2= a = omega^2 sep => omega^2 = G mu / sep^3 => v_i = r_i omega = r_i sqrt(G (m1+ m2) / sep^3); [r_i G (m1+m2) / sep^3]) = kpc (km^2 kpc Msun / (Msun s^2 kpc^3))^1/2 = km/s
    // [omega] = [sqrt(G (m1+ m2) / sep^3)] = [km^2 kpc Msun / Msun s^2 kpc^3]^1/2 = km / s kpc
    let omega = (GG*(m1+m2)/(seperation*seperation*seperation)).sqrt();//.3; //1.1 fudge factor to make it a tad elliptical
    println!("ensure timestep less than 2/omega = {}", (2.0*KM_IN_KPC/omega));

    let bh_1 = BlackHole {
        mass: m1,
        position: [offset, 0.0, 0.0],
        velocity: [0.0, offset * omega, 0.0],
        acceleration: [0.0; 3],
        jerk: [0.0; 3],
        snap: [0.0; 3],
    };
    let bh_2 = BlackHole {
        mass: m2,
        position: [offset - seperation, 0.0, 0.0],
        velocity: [0.0, (offset - seperation) * omega, 0.0],
        acceleration: [0.0; 3],
        jerk: [0.0; 3],
        snap: [0.0; 3],
    };
    vec![bh_1, bh_2]
}

fn compute_velocity_dispersion(v_points: &Vec<f64>, p_of_v_given_r: &Vec<Vec<f64>>) -> Vec<f64> {
    // note for isotropic systems <v_i> = 0 so sigma_i = <v_i^2> - <v_i>^2 = <v_i>^2 so then <v^2> = sum <v_i^2> = sum sigma_i^2 = sigma_tot^2
    let mut v_disp: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];

    for i in 0..SPACIAL_GRID_NUM {
        for j in 0..VELOCITY_GRID_NUM {
            v_disp[i] += v_points[j] * v_points[j] * p_of_v_given_r[i][j];
        }
        // 1D vdisp is 1/3 rd of 3D vdisp
        v_disp[i] /= 3.0
    }

    v_disp
}

fn recover_rho(f: &Vec<Vec<f64>>, v_points: &Vec<f64>) -> Vec<f64> {
    let mut rho: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];

    for i in 0..SPACIAL_GRID_NUM {
        rho[i] = f[i][0] * 4.0 * PI * pow(0.5*v_points[1], 3)/3.0;
        for j in 1..(VELOCITY_GRID_NUM - 1) {
            rho[i] += f[i][j] * 4.0 * PI * (pow(0.5*(v_points[j+1] + v_points[j]), 3) - pow(0.5*(v_points[j] + v_points[j-1]), 3)) / 3.0;
        }
        // Assumes linear spacing or r for last step, shouldn't matter since rho should basically be zero out here anyway
        rho[i] += f[i][VELOCITY_GRID_NUM - 1] * 4.0 * PI * (pow(1.5*v_points[VELOCITY_GRID_NUM - 1] - 0.5*v_points[VELOCITY_GRID_NUM - 2], 3) - pow(0.5*(v_points[VELOCITY_GRID_NUM - 1] + v_points[VELOCITY_GRID_NUM-2]), 3)) / 3.0;
    }

    rho
}

fn plot_rho_recovery_test(r_points: &Vec<f64>, rho_points: &Vec<f64>, recovered_rho_points: &Vec<f64>, filename: &str) -> Result<(), Box<dyn std::error::Error>> { 
    let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;

    for i in 0..rho_points.len() {
        if rho_points[i] > y_max {
            y_max = rho_points[i]
        }
        if rho_points[i] < y_min {
            y_min = rho_points[i]
        }
    }

    let mut chart = ChartBuilder::on(&root)
        .caption("Invertability test for eddington inversion", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d((r_points[0]..r_points[r_points.len() - 1]),
            ((y_min+1e-4) * match y_min.signum() {
                1.0 => 0.9,
                -1.0 => 1.1,
                _ => panic!("number has no sign, is probably NaN")
            }..y_max * match y_max.signum() {
                1.0 => 1.1,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN")
            })
        )?;

    chart.configure_mesh()
        .x_desc("r (kpc)")           // X-axis label
        .y_desc("rho (M_sun / kpc^3)") // Y-axis label
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    let mut plot_profile: Vec<(f64, f64)> = (0..r_points.len())
        .map(|i| (r_points[i], rho_points[i]))
        .collect();

    chart.draw_series(LineSeries::new(plot_profile, &BLUE))?.label("input density profile").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    let mut plot_profile: Vec<(f64, f64)> = (0..r_points.len())
        .map(|i| (r_points[i], recovered_rho_points[i]))
        .collect();

    chart.draw_series(LineSeries::new(plot_profile, &RED))?.label("density profile calculated by eddington inversion").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    root.present()?;
    println!("rho recovery plot saved as {}", filename);
    Ok(())
}