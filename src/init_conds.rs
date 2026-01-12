use core::f64;
use std::f64::consts::PI;
use statrs::function::gamma as Gamma;

use crate::particle::Particle;
use crate::forces::GG;
use crate::eddington_inverter::*;
use crate::plotting::plot_check_function;
use crate::plotting::plot_function;
use crate::vectors::*;
use crate::logging::load_file;

pub const SPACIAL_GRID_NUM: usize = 10000;
pub const VELOCITY_GRID_NUM: usize = 10000;
const MASS_ERROR_TOLERANCE: f64 = 0.001;

pub fn combine(file_paths: Vec<String>, offsets: Vec<[f64; 3]>, velocity_offsets: Vec<[f64; 3]>) -> Vec<Particle> {
    if file_paths.len() != offsets.len() {
        panic!("The number of constituent simulations ({}) must match the number of offsets ({})", file_paths.len(), offsets.len())
    }

    if file_paths.len() != velocity_offsets.len() {
        panic!("The number of constituent simulations ({}) must match the number of velocity offsets ({})", file_paths.len(), velocity_offsets.len())
    }

    let mut combo_data: Vec<Particle> = Vec::new();

    for n in 0..file_paths.len() {
        let mut ingredient = load_file(file_paths[n].clone()).unwrap().data;
        for i in 0..ingredient.len() {
            ingredient[i].position = add(&ingredient[i].position, &offsets[n]);
            ingredient[i].velocity = add(&ingredient[i].velocity, &velocity_offsets[n]);
        }
        combo_data.append(&mut ingredient);
    }

    combo_data
}

pub fn abg_profile_init_conds(alpha: f64, beta: f64, gamma: f64, r_s: f64, rho_s: f64, r_cutoff_option: Option<f64>, particle_num: usize, output_path: String) -> Vec<Particle> {

    if gamma > 3.0 {
        panic!("Total mass diverges at small r");
    }

    let cuspy: bool = gamma > 0.0;

    let particles: Vec<Particle> = match r_cutoff_option {
        Some(r_cutoff) => {
            let rho_cutoff = rho_s / ((r_cutoff / r_s).powf(gamma) * (1.0 + (r_cutoff / r_s).powf(alpha)).powf((beta - gamma) / alpha));
            let r_decay = 0.3 * r_cutoff;
            let delta = (r_cutoff / r_decay) - ((gamma + beta * (r_cutoff / r_s).powf(alpha)) / (1.0 + (r_cutoff / r_s).powf(alpha)));

            let i_m = inner_mass_integral(r_cutoff/ r_s, alpha, beta, gamma);
            let i_m_cut = outer_mass_integral(r_s, alpha, beta, gamma, r_cutoff, delta, r_decay, 10.0 * r_decay);

            let total_mass = 4.0 * PI * r_s.powi(3) * rho_s * (i_m + i_m_cut);
            if total_mass.is_infinite() || total_mass.is_nan() {
                panic!("Total mass diverges, check your parameters!");
            }
            if total_mass < 0.0 {
                panic!("Total mass negative, check your parameters!");
            }

            let particle_mass = total_mass / (particle_num as f64);

            let r_max = {
                let mut r_max = 0.0;
                let mut enclosed_mass_frac = 0.0;
                let mut i_m;
                while enclosed_mass_frac < 1.0 - MASS_ERROR_TOLERANCE {
                    r_max += r_s;
                    match r_max < r_cutoff {
                        true => i_m = inner_mass_integral(r_max/ r_s, alpha, beta, gamma),
                        false => i_m = inner_mass_integral(r_cutoff/ r_s, alpha, beta, gamma) + outer_mass_integral(r_s, alpha, beta, gamma, r_cutoff, delta, r_decay, r_max)
                    }
                    enclosed_mass_frac = (4.0 * PI * r_s.powi(3) * rho_s * i_m) / total_mass;
                }
                r_max
            };

            let rho = |r: f64| -> f64 {
                if r < r_cutoff {
                    rho_s / ((r / r_s).powf(gamma) * (1.0 + (r / r_s).powf(alpha)).powf((beta - gamma) / alpha))
                } else {
                    rho_cutoff * (r / r_cutoff).powf(delta) * ((r_cutoff - r) / r_decay).exp()
                }
            };

            let tail = AsymptoticTail::Exponential(GG * total_mass / r_decay, -delta - 3.0);
            // [GG * M / r_decay] = [km^2 kpc / Msun s^2] [M_sun] / [kpc] = [km^2 / s^2]

            let (r_points, p_of_r, v_points, p_of_v_given_r) = generate_init_conds_from_rho(rho, &r_max, &cuspy, tail, output_path.clone());

            let particles = sample_from_distribution(&r_points, &p_of_r, &v_points, &p_of_v_given_r, particle_num, particle_mass);
    
            check_output(&particles, &rho, output_path);

            particles
        },
        None => {
            if beta <= 3.0 {
                panic!("Total mass diverges at large r, please provide a cutoff radius");
            }
            let i_m = Gamma::gamma((beta - 3.0)/ alpha) * Gamma::gamma((3.0 - gamma)/ alpha) / (alpha * Gamma::gamma((beta - gamma)/ alpha));
            let total_mass = 4.0 * PI * r_s.powi(3) * rho_s * i_m ;
            if total_mass.is_infinite() || total_mass.is_nan() {
                panic!("Total mass diverges, check your parameters!");
            }
            if total_mass < 0.0 {
                panic!("Total mass negative, check your parameters!");
            }

            let particle_mass = total_mass / (particle_num as f64);

            let r_max = {
                let mut r_max = 0.0;
                let mut enclosed_mass_frac = 0.0;
                while enclosed_mass_frac < 1.0 - MASS_ERROR_TOLERANCE {
                    r_max += r_s;
                    enclosed_mass_frac = (4.0 * PI * r_s.powi(3) * rho_s * inner_mass_integral(r_max/ r_s, alpha, beta, gamma)) / total_mass;
                }
                r_max
            };

            let rho = |r: f64| -> f64 {
                rho_s / ((r / r_s).powf(gamma) * (1.0 + (r / r_s).powf(alpha)).powf((beta - gamma) / alpha))
            };

            let tail = AsymptoticTail::PowerLaw(beta);

            let (r_points, p_of_r, v_points, p_of_v_given_r) = generate_init_conds_from_rho(rho, &r_max, &cuspy, tail, output_path.clone());

            let particles = sample_from_distribution(&r_points, &p_of_r, &v_points, &p_of_v_given_r, particle_num, particle_mass);
    
            check_output(&particles, &rho, output_path);

            particles
        }
    };

    particles
}

fn inner_mass_integral(q: f64, alpha: f64, beta: f64, gamma: f64) -> f64 {
    // Numerical integration for general case
    // But for specific αβγ combinations, we have analytic forms
    
    if alpha == 1.0 && beta == 3.0 && gamma == 1.0 {
        // NFW: ∫_0^q x/(1+x)² dx
        return (1.0 + q).ln() - q / (1.0 + q);
    } else if alpha == 2.0 && beta == 5.0 && gamma == 0.0 {
        // Plummer: ∫_0^q x²/(1+x²)^{5/2} dx
        return q.powi(3) / (3.0 * (q.powi(2) + 1.0).powf(1.5));
    } else {
        // General case: numerical integration
        let n_points = 1000;
        let mut integral = 0.0;
        let dx = q / (n_points as f64);
        for i in 0..(n_points - 1) {
            let x = q * (i as f64) / (n_points as f64);
            let x_2 = q * ((i+1) as f64) / (n_points as f64);
            
            let integrand = x.powf(2.0 - gamma) / 
                           (1.0 + x.powf(alpha)).powf((beta - gamma) / alpha);
            let integrand_2 = x_2.powf(2.0 - gamma) / 
                           (1.0 + x_2.powf(alpha)).powf((beta - gamma) / alpha);
            
                           integral += 0.5 * (integrand + integrand_2) * dx;
        }
        return integral;
    }
}

fn outer_mass_integral(r_s: f64, alpha: f64, beta: f64, gamma: f64, r_cutoff: f64, delta: f64, r_decay: f64, r_max: f64) -> f64 {
    let mut integral = 0.0;
    let t_max = r_max - r_cutoff;
    let n_points = 1000;
    let dt = t_max / (n_points as f64);
    for i in 0..n_points {
        let t = t_max * (i as f64) / (n_points as f64);
        let t_2 = t_max * ((i+1) as f64) / (n_points as f64);
        
        let r = r_cutoff + t;
        let r_2 = r_cutoff + t_2;
        
        let integrand = r.powi(2) * (r / r_cutoff).powf(delta) * 
                        (-t / r_decay).exp();
        let integrand_2 = r_2.powi(2) * (r_2 / r_cutoff).powf(delta) * 
                        (-t_2 / r_decay).exp();

        integral += 0.5 * (integrand + integrand_2) * dt;
    }
    let q = r_cutoff / r_s;
    integral /= r_s.powi(3) * q.powf(gamma) * (1.0 + q.powf(alpha)).powf((beta - gamma)/alpha);
    integral
}

pub fn plummer_init_conds(r_s: f64, total_mass: f64, particle_num: usize, output_path: String) -> Vec<Particle> {
    
    let rho_s = 3.0 * total_mass / (4.0 * PI * r_s.powi(2));

    let particles = abg_profile_init_conds(2.0, 5.0, 0.0, r_s, rho_s, None, particle_num, output_path);

    particles
}

pub fn nfw_init_conds(r_s: f64, rho_s: f64, r_cutoff: f64, particle_num: usize, output_path: String) -> Vec<Particle> {
    
    let particles = abg_profile_init_conds(1.0, 3.0, 1.0, r_s, rho_s, Some(r_cutoff), particle_num, output_path);

    particles
}

fn generate_init_conds_from_rho<T: Fn(f64) -> f64>(rho: T, r_max: &f64, cuspy: &bool, tail: AsymptoticTail, output_path: String) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
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

    let rho_points = generate_rho_points(&rho, &r_points);

    let potential_points = generate_potential_points(&rho_points, &r_points);
    let v_max = (-2.0 * potential_points[0]).sqrt();
    
    //plot_check_function(&r_points, &|r: f64| -> f64 {- GG*1e8/((1.0 + r.powi(2)).sqrt())}, &potential_points, &"Potential_Points_test.png", &"Potential Points Check", &"r (kpc)", &"V (M_sun km^2/s^2)")?;

    let mut v_points: Vec<f64> = vec![0.0; VELOCITY_GRID_NUM];
    for j in 0..VELOCITY_GRID_NUM {
        v_points[j] = (j as f64) * v_max / ((VELOCITY_GRID_NUM - 1) as f64);
    }

    let mut p_of_r: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    let mut total_mass: f64 = 0.0;

    p_of_r[0] = rho_points[0] * 4.0 * PI * (0.5*r_points[1]).powi(3)/3.0;
    total_mass += p_of_r[0];
    for i in 1..(SPACIAL_GRID_NUM - 1) {
        p_of_r[i] = rho_points[i] * 4.0 * PI * ((0.5*(r_points[i+1] + r_points[i])).powi(3) - (0.5*(r_points[i] + r_points[i-1])).powi(3)) / 3.0;
        total_mass += p_of_r[i];
    }
    // Assumes linear spacing of r for last step, shouldn't matter since rho should basically be zero out here anyway
    p_of_r[SPACIAL_GRID_NUM - 1] = rho_points[SPACIAL_GRID_NUM - 1] * 4.0 * PI * ((1.5*r_points[SPACIAL_GRID_NUM - 1] - 0.5*r_points[SPACIAL_GRID_NUM - 2]).powi(3) - (0.5*(r_points[SPACIAL_GRID_NUM - 1] + r_points[SPACIAL_GRID_NUM-2])).powi(3)) / 3.0;
    total_mass += p_of_r[SPACIAL_GRID_NUM - 1];

    //Normalizing
    for i in 0..SPACIAL_GRID_NUM {
        p_of_r[i] /= total_mass;
    }

    let f: Vec<Vec<f64>> = compute_phase_space_density(&rho_points, &v_points, &potential_points, &cuspy, tail);

    let mut p_of_v_given_r: Vec<Vec<f64>> = vec![vec![0.0; VELOCITY_GRID_NUM]; SPACIAL_GRID_NUM];
    let mut normalization_check: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    
    let mut normalization_min: f64 = f64::MAX;
    let mut normalization_max: f64 = 0.0;
    for i in 0..SPACIAL_GRID_NUM {
        let v_escape = (-2.0 * potential_points[i]).sqrt();
        let v_esc_arg = ((v_escape/v_max) * VELOCITY_GRID_NUM as f64) as usize - 1;

        p_of_v_given_r[i][0] = (4.0 * PI / 3.0) * (f[i][0] / rho_points[i]) * (0.5*v_points[1]).powi(3);
        normalization_check[i] += p_of_v_given_r[i][0];
        
        for j in 1..v_esc_arg {
            p_of_v_given_r[i][j] = (4.0 * PI / 3.0) * (f[i][j] / rho_points[i]) * ((0.5*(v_points[j+1] + v_points[j])).powi(3) - (0.5*(v_points[j] + v_points[j-1])).powi(3));
            normalization_check[i] += p_of_v_given_r[i][j];
        }

        normalization_min = normalization_min.min(normalization_check[i]);
        normalization_max = normalization_max.max(normalization_check[i]);

        //renormalizing
        for j in 0..VELOCITY_GRID_NUM {
            p_of_v_given_r[i][j] /= normalization_check[i];
        }
    }

    if normalization_min < 0.9 || normalization_max > 1.1 {            
            eprintln!("Velocity space distribution is not normalized, this likely indicates an issue with the asymptotic behavior at large r. Please consult the diagnostic plots to see if the accuracy is sufficient for your needs. Normalization ∈ [{normalization_min}, {normalization_max}]")
        }
    
    plot_function(&r_points, &normalization_check, (output_path.clone() + &"/Normalization.png").as_str(), &"Velocity Normalization Check", &"Normalization", &"r").unwrap();
    plot_function(&r_points, &compute_velocity_dispersion(&v_points, &p_of_v_given_r), (output_path.clone() + &"/Eddington_Inversion_Velocity_Dispersion_test.png").as_str(), &"Velocity Dispersion Reconstructed From Eddington Inversion", &"r (kpc)", &"sigma^2 (km^2/s^2)").unwrap();
    plot_check_function(&r_points, &rho, &recover_rho_from_f(&f, &v_points), (output_path.clone() + &"/Eddington_Inversion_Density_test.png").as_str(), &"Density Reconstructed From Eddington Inversion", &"r (kpc)", &"rho (M_sun / kpc^3)").unwrap();
    
    (r_points, p_of_r, v_points, p_of_v_given_r)
}

fn sample_from_distribution(r_points: &Vec<f64>, p_of_r: &Vec<f64>, v_points: &Vec<f64>, p_of_v_given_r: &Vec<Vec<f64>>, particle_num: usize, particle_mass: f64) -> Vec<Particle> {
    let mut particles: Vec<Particle> = Vec::with_capacity(particle_num);

    let cumulative_p_of_r: Vec<f64> = {
        let mut cumulative: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
        cumulative[0] = p_of_r[0];
        for i in 1..SPACIAL_GRID_NUM {
            cumulative[i] = cumulative[i-1] + p_of_r[i];
        }
        cumulative
    };

    let cumulative_p_of_v_given_r: Vec<Vec<f64>> = {
        let mut cumulative: Vec<Vec<f64>> = vec![vec![0.0; VELOCITY_GRID_NUM]; SPACIAL_GRID_NUM];
        for i in 0..SPACIAL_GRID_NUM {
            cumulative[i][0] = p_of_v_given_r[i][0];
            for j in 1..VELOCITY_GRID_NUM {
                cumulative[i][j] = cumulative[i][j-1] + p_of_v_given_r[i][j];
            }
        }
        cumulative
    };

    for _ in 0..particle_num {
        let r_rand = rand::random::<f64>();
        let r_index: usize;
        
        let v_rand = rand::random::<f64>();

        let r = {
            let mut high = SPACIAL_GRID_NUM - 1;
            let mut low = 0;
            while high - low > 1 {
                let mid = (low + high) / 2;

                //check if above or below
                match r_rand > cumulative_p_of_r[mid] {
                    true => low = mid,
                    false => high = mid
                };
            }
            let t = (r_rand - cumulative_p_of_r[low]) / (cumulative_p_of_r[high] - cumulative_p_of_r[low]);
            r_index = low;
            r_points[low] * (1.0 - t) + r_points[high] * t
        };

        let v = {
            let mut high = VELOCITY_GRID_NUM - 1;
            let mut low = 0;
            while high - low > 1 {
                let mid = (low + high) / 2;

                //check if above or below
                match v_rand > cumulative_p_of_v_given_r[r_index][mid] {
                    true => low = mid,
                    false => high = mid
                };
            }
            let t = (v_rand - cumulative_p_of_v_given_r[r_index][low]) / (cumulative_p_of_v_given_r[r_index][high] - cumulative_p_of_v_given_r[r_index][low]);
            v_points[low] * (1.0 - t) + v_points[high] * t
        };

        particles.push(Particle {
            mass: particle_mass,
            position: scalar_multiply(&r, &random_unit_vector()),
            velocity: scalar_multiply(&v, &random_unit_vector()),
            acceleration: [0.0; 3],
            jerk: [0.0; 3],
            snap: [0.0; 3],
        });
    }

    particles
}

fn random_unit_vector() -> [f64; 3] {
    let phi = rand::random::<f64>() * 2.0 * PI; // phi uniform in [0, 2pi]
    let cos_theta = (rand::random::<f64>() * 2.0) - 1.0; // cos(theta) uniform in [-1, 1]
    let theta = cos_theta.acos();

    let x = theta.sin() * phi.cos();
    let y = theta.sin() * phi.sin();
    let z = theta.cos();

    [x, y, z]
}

pub fn binary_init_conds(m1: f64, m2: f64, min_seperation: f64, eccentricity: f64, initial_seperation: f64) -> Vec<Particle> {
    // semi-major axis
    let a = min_seperation / (1.0 - eccentricity);
    
    // Check validity
    if initial_seperation < min_seperation {
        panic!("Your initial seperation may not be less than your minimum seperation")
    }
    if a > 0.0 && initial_seperation > 2.0*a {
        panic!("Your initial seperation exceeds the maximum seperation of your elliptical orbit")
    }

    // Offset so centered on COM
    // 0 = m1*x + m2 (x - sep) = - m2 sep + (m1 + m2)x => x = m2 sep/ (m1 + m2)
    let offset = m2 * initial_seperation / (m1+ m2);
    
    // v_rel = sqrt(G (m1 + m2) (2/r - 1/a))
    let relative_velocity = (GG * (m1 + m2) * (2.0/initial_seperation - 1.0/a)).sqrt();

    // m1 v1 + m2 v2 = 0 & v1 - v2 = v_rel => m1 v1 + m2 (v1 - v_rel) = 0 => v1 = m2 v_rel / (m1 + m2)
    // => v2 = - m1 v_rel / (m1 + m2)
    let v1 = m2 * relative_velocity / (m1 + m2);
    let v2 = - m1 * relative_velocity / (m1 + m2);
    let momentum = m1 * m2 * relative_velocity / (m1 + m2);
    
    // In COM frame:
    // L = r1 x p1 + r2 x p2 = r1 x p1 - r2 x p1 = (r1 - r2) x p1 = sep x p
    // => L = min_sep * pmax
    let vmax = (GG * (m1 + m2) * (2.0/min_seperation - 1.0/a)).sqrt();
    let pmax = m1 * m2 * vmax / (m1 + m2);
    let angular_momentum = min_seperation * pmax;

    // => L = sep * p * sin(phi)
    // => phi = arcsin(L/(sep*p))
    let phi = (angular_momentum / (initial_seperation * momentum)).asin();

    let particle_1 = Particle {
        mass: m1,
        position: [offset, 0.0, 0.0],
        velocity: [-v1 * phi.cos(), v1 * phi.sin(), 0.0],
        acceleration: [0.0; 3],
        jerk: [0.0; 3],
        snap: [0.0; 3],
    };
    let particle_2 = Particle {
        mass: m2,
        position: [offset - initial_seperation, 0.0, 0.0],
        velocity: [-v2 * phi.cos(), v2 * phi.sin(), 0.0],
        acceleration: [0.0; 3],
        jerk: [0.0; 3],
        snap: [0.0; 3],
    };
    vec![particle_1, particle_2]
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

fn recover_rho_from_f(f: &Vec<Vec<f64>>, v_points: &Vec<f64>) -> Vec<f64> {
    let mut rho: Vec<f64> = vec![0.0; SPACIAL_GRID_NUM];
    
    for i in 0..SPACIAL_GRID_NUM {
        rho[i] = f[i][0] * 4.0 * PI * (0.5*v_points[1]).powi(3)/3.0;
        for j in 1..(VELOCITY_GRID_NUM - 1) {
            rho[i] += f[i][j] * 4.0 * PI * ((0.5*(v_points[j+1] + v_points[j])).powi(3) - (0.5*(v_points[j] + v_points[j-1])).powi(3)) / 3.0;
        }
        // Assumes linear spacing or r for last step, shouldn't matter since rho should basically be zero out here anyway
        rho[i] += f[i][VELOCITY_GRID_NUM - 1] * 4.0 * PI * ((1.5*v_points[VELOCITY_GRID_NUM - 1] - 0.5*v_points[VELOCITY_GRID_NUM - 2]).powi(3) - (0.5*(v_points[VELOCITY_GRID_NUM - 1] + v_points[VELOCITY_GRID_NUM-2])).powi(3)) / 3.0;
    }
    rho
}

fn check_output<T: Fn(f64) -> f64>(particles: &Vec<Particle>, rho_analytic: &T, output_path: String) {
    let bin_count = (particles.len() / 1000).max(5);
    let mut rho: Vec<f64> = vec![0.0; bin_count];
    let mut v_disp: Vec<f64> = vec![0.0; bin_count];
    let mut counts: Vec<usize> = vec![0; bin_count];

    let (r_min, r_max): (f64, f64) = {
        let mut r_max = 0.0;
        let mut r_min = f64::MAX;
        for p in particles.iter() {
            let r = magnitude(&p.position);
            if r > r_max {
                r_max = r;
            }
            if r < r_min && r > 0.0 {
                r_min = r;
            }
        }
        (r_min*0.95, r_max*1.05)
    };

    let log_bin_width = (r_max.ln() - r_min.ln()) / (bin_count as f64);

    let bin_edges: Vec<f64> = {
        let mut r_points: Vec<f64> = vec![0.0; bin_count + 1];
        for i in 0..bin_count + 1 {
            r_points[i] = (r_min.ln() + (i as f64) * (r_max.ln() - r_min.ln()) / ((bin_count) as f64)).exp()
        }
        r_points
    };
    let (r_points, bin_volumes): (Vec<f64>, Vec<f64>) = {
        let mut r_points: Vec<f64> = vec![0.0; bin_count];
        let mut bin_volumes: Vec<f64> = vec![0.0; bin_count];
        for i in 0..bin_count {
            r_points[i] = (bin_edges[i] * bin_edges[i+1]).sqrt(); // geometric mean for log spaced bins
            bin_volumes[i] = (4.0/3.0) * PI * (bin_edges[i+1].powi(3) - bin_edges[i].powi(3));
        }
        (r_points, bin_volumes)
    };
    
    for p in particles.iter() {
        let r = magnitude(&p.position);
        if r > r_max || r < r_min {
            continue;
        }
        
        // Find bin index directly (no binary search needed for uniform log bins)
        let bin_index = ((r.ln() - r_min.ln()) / log_bin_width).floor() as usize;
        rho[bin_index] += p.mass / bin_volumes[bin_index];
        let v2 = dot_product(&p.velocity, &p.velocity);
        v_disp[bin_index] += v2;
        counts[bin_index] += 1;
    }

    for i in 0..rho.len() {
        if counts[i] > 0 {
            v_disp[i] /= counts[i] as f64;
            v_disp[i] /= 3.0; // 1D velocity dispersion
        }
    }

    plot_check_function(&r_points, rho_analytic, &rho, &(output_path.clone() + &"/Sampler_Density_test.png"), &"Sampler Density Test", &"r (kpc)", &"rho (M_sun / kpc^3)").unwrap();
    plot_function(&r_points , &v_disp, &(output_path.clone() + &"/Sampler_Velocity_Dispersion_test.png"), &"Sampler Velocity Dispersion Test", &"r (kpc)", &"v_disp (km/s)").unwrap();

}