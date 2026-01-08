use plotters::prelude::*;
use statrs::statistics::Statistics;
use core::f64;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use crate::particle::Particle;
use crate::vectors::magnitude;
use crate::forces::{calculate_energy, calculate_momentum, calculate_angular_momentum};
use crate::logging::load_checkpoint;

pub fn plot_trajectories(data_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let scale = 4;
    let binding = &(filename.to_owned() + ".png");
    let root = BitMapBackend::new(binding, (1024*scale, 1024*scale)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut trajectories: Vec<Vec<(f64, f64, f64)>> = Vec::new();
    let mut time: Vec<f64> = Vec::new();
    let mut bound: f64 = 0.0;
    
    let mut file_num: usize = 0;
    let mut last_data: Vec<Particle> = Vec::new();

    while Path::new(&format!("{}/snapshot_{:03}.log", data_directory, file_num)).exists() {
        let state = load_checkpoint(data_directory, &file_num)?;
        let data = state.data;
        time.push(state.time);

        if trajectories.is_empty() {
            trajectories = vec![Vec::new(); data.len()];
        }

        let mut ave_dist = 0.0;
        // Add current positions to each particle's trajectory
        for (particle_idx, particle) in data.iter().enumerate() {
            trajectories[particle_idx].push((
                particle.position[0],
                particle.position[1],
                particle.position[2],
            ));

            ave_dist += magnitude(&particle.position);
        }
        ave_dist /= data.len() as f64;

        bound = bound.max(ave_dist);

        
        last_data = data;
        file_num += 1;
    }

    // Movie
    let mut movie_file_path = PathBuf::from(filename);
    let file_name = movie_file_path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| "default".to_string());
    movie_file_path.pop();
    movie_file_path.push("movie_dump");
    fs::create_dir_all(&movie_file_path)?;
    movie_file_path.push(file_name);
    let movie_file_name = movie_file_path.to_string_lossy().into_owned();
    
    for i in 0..trajectories.len() {
        let binding = &(movie_file_name.clone() + &format!("_{:03}.png", i));
        let instantaneous_positions: Vec<(f64, f64, f64)> = trajectories.iter().map(|trajectory| {trajectory[i]}).collect();
        plot_positions(&instantaneous_positions, bound, scale, binding, time[i])?
    }

    if Command::new("ffmpeg").arg("-version").output().is_err() {
        return Err("ffmpeg is not installed or not found in PATH".into())
    }

    let movie_name = &(filename.to_owned() + ".mkv");

    if Path::new(movie_name).exists() {
        fs::remove_file(movie_name)?;
    }

    Command::new("ffmpeg")
        .args(&[
            "-framerate",  "12",
            "-i", &(movie_file_name + "_%03d.png"),
            "-q:v", "0",
            movie_name,
        ]).status()?;

    // Static trails plot

    let mut chart = ChartBuilder::on(&root)
        .caption("Trajectories", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_3d(- bound..bound, - bound..bound, - bound..bound)?;

    chart.configure_axes().draw()?;

    let colors = vec![&RED, &BLUE, &GREEN, &YELLOW, &CYAN, &MAGENTA];

    let particle_count = trajectories.len();

    // Plot trajectories for each particle
    for (particle_idx, traj) in trajectories.iter().enumerate() {
        if traj.is_empty() {
            continue;
        }
        
        let color = colors[particle_idx % colors.len()];
        
        chart.draw_series(LineSeries::new(
            traj.clone(),
            color.stroke_width((20 / particle_count.isqrt() as u32).max(1)),
        ))?
        //.label(format!("Particle {}", particle_idx + 1))
        //.legend(move |(x, y)| {PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))})
        ;
    }

    // Add final positions as markers
    for (particle_idx, particle) in last_data.iter().enumerate() {
        let color = colors[particle_idx % colors.len()];
        
        chart.draw_series(PointSeries::of_element(
            vec![(particle.position[0], particle.position[1], particle.position[2])],
            (50 / particle_count.isqrt() as i32).max(2),
            color.filled(),
            &|c, s, st| {
                return EmptyElement::at(c) + Circle::new((0, 0), s, st);
            },
        ))?;
    }

    //chart.configure_series_labels().background_style(&WHITE.mix(0.8)).border_style(&BLACK).draw()?;

    root.present()?;
    println!("Plot saved as {}", filename);
    Ok(())
}

fn plot_positions(positions: &Vec<(f64, f64, f64)>, bound: f64, scale: u32, filename: &str, time: f64) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (1024*scale, 1024*scale)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Trajectories t={:03}", time), ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_3d(- bound..bound, - bound..bound, - bound..bound)?;

    chart.configure_axes().draw()?;

    let colors = vec![&RED, &BLUE, &GREEN, &YELLOW, &CYAN, &MAGENTA];

    let particle_count = positions.len();

    for (particle_idx, position) in positions.iter().enumerate() {
        let color = colors[particle_idx % colors.len()];
        
        chart.draw_series(PointSeries::of_element(
            vec![*position],
            (50 / particle_count.isqrt() as i32).max(2),
            color.filled(),
            &|c, s, st| {
                return EmptyElement::at(c) + Circle::new((0, 0), s, st);
            },
        ))?;
    }

    root.present()?;
    Ok(())
}

pub fn plot_energy(data_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (2*1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((1,2));

    let mut kinetic: Vec<(f64, f64)> = Vec::new();
    let mut potential: Vec<(f64, f64)> = Vec::new();
    let mut tot_energy: Vec<(f64, f64)> = Vec::new();
    
    let mut file_num = 0; 
    let mut max_time = 0.0;
    let mut y_max = f64::MIN;
    let mut y_min = f64::MAX;
    while Path::new(&format!("{}/snapshot_{:03}.log", data_directory, file_num)).exists() {
        let state = load_checkpoint(data_directory, &file_num)?;

        let (kinetic_energy, potential_energy) = calculate_energy(&state.data);

        kinetic.push((state.time, kinetic_energy));
        potential.push((state.time, potential_energy));
        tot_energy.push((state.time, kinetic_energy + potential_energy));

        y_max = y_max.max(kinetic_energy);
        y_min = y_min.min(potential_energy);

        max_time = state.time;
        file_num += 1
    }
    y_max = y_max.max(1e-3);
    y_min = y_min.min(-1e-3);
    
    let mut metric_factor = 0;
    while max_time < (10.0 as f64).powi(metric_factor) {
        metric_factor -= 3
    }
    for i in 0..kinetic.len() {
        kinetic[i].0 /= (10.0 as f64).powi(metric_factor);
        potential[i].0 /= (10.0 as f64).powi(metric_factor);
        tot_energy[i].0 /= (10.0 as f64).powi(metric_factor);
    }
    let time_unit = match metric_factor {
        0 => "Gyr",
        -3 => "Myr",
        -6 => "Kyr",
        -9 => "Year",
        _ => &format!("1e{} Years", 9 - metric_factor) // largest supported unit
    };

    let (conservation_check, conservation_bounds) = {
        let mut conservation_check = Vec::with_capacity(tot_energy.len());
        let mut y_min = f64::MAX;
        let mut y_max = f64::MIN;
        for i in 0..tot_energy.len() {
            let non_conservation_amount = (tot_energy[i].1 - tot_energy[0].1) / tot_energy[0].1.abs().max(1e-10); // avoid division by 0
            conservation_check.push((tot_energy[i].0, non_conservation_amount));
            y_min = y_min.min(non_conservation_amount);
            y_max = y_max.max(non_conservation_amount);
        }
        y_max = y_max.max(1e-15);
        y_min = y_min.min(-1e-15);
        (conservation_check, (y_min, y_max))
    };

    let mut energy_chart = ChartBuilder::on(&areas[0])
        .caption("Total Energy", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..max_time, y_min * 1.1..y_max * 1.1)?;

    energy_chart.configure_mesh()
        .x_desc(format!("Time ({time_unit})"))           
        .y_desc("Energy (M_sun km^2 / s^2)")
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    energy_chart.draw_series(LineSeries::new(kinetic, &BLUE))?.label("Kinetic").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    energy_chart.draw_series(LineSeries::new(potential, &GREEN))?.label("Potential").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
    energy_chart.draw_series(LineSeries::new(tot_energy, &BLACK))?.label("Total").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    energy_chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    let mut conservation_chart = ChartBuilder::on(&areas[1])
        .caption("Energy Conservation", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..max_time, conservation_bounds.0 * 1.1..conservation_bounds.1 * 1.1)?;

    conservation_chart.configure_mesh()
        .x_desc(format!("Time ({time_unit})"))           
        .y_desc("(E - E_0) / |E_0|")
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    conservation_chart.draw_series(LineSeries::new(conservation_check, &BLACK))?;

    root.present()?;
    println!("Energy plot saved as {}", filename);
    Ok(())
}

pub fn plot_momentum(data_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (2*1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((1,2));

    let mut momentum: Vec<Vec<(f64, f64)>> = vec![Vec::new(); 3];
    
    let mut file_num = 0; 
    let mut max_time = 0.0;
    let mut y_max = f64::MIN;
    let mut y_min = f64::MAX;
    while Path::new(&format!("{}/snapshot_{:03}.log", data_directory, file_num)).exists() {
        let state = load_checkpoint(data_directory, &file_num)?;

        let momentum_data = calculate_momentum(&state.data);

        for i in 0..3 {
            momentum[i].push((state.time, momentum_data[i]));
        }
        
        y_max = y_max.max(momentum_data.max());
        y_min = y_min.min(momentum_data.min());

        max_time = state.time;
        file_num += 1
    }
    y_max = y_max.max(1e-3);
    y_min = y_min.min(-1e-3);

    let mut metric_factor = 0;
    while max_time < (10.0 as f64).powi(metric_factor) {
        metric_factor -= 3
    }
    for n in 0..momentum[0].len() {
        for i in 0..3 {
            momentum[i][n].0 /= (10.0 as f64).powi(metric_factor);
        }
    }
    let time_unit = match metric_factor {
        0 => "Gyr",
        -3 => "Myr",
        -6 => "Kyr",
        -9 => "Year",
        _ => &format!("1e{} Years", 9 - metric_factor) // largest supported unit
    };

    let (conservation_check , conservation_bounds) = {
        let mut conservation_check = vec![Vec::with_capacity(momentum[0].len()), Vec::with_capacity(momentum[0].len()),Vec::with_capacity(momentum[0].len())];
        let mut y_min = f64::MAX;
        let mut y_max = f64::MIN;
        for n in 0..momentum[0].len() {
            let non_conservation_amount = {
                let mut non_conservation_amount = Vec::new();
                for i in 0..3 {
                    non_conservation_amount.push((momentum[i][n].1 - momentum[i][0].1) / momentum[i][0].1.abs().max(1e-10)) // to avoid division by zero
                }
                non_conservation_amount
            };

            y_min = y_min.min(non_conservation_amount.clone().min());
            y_max = y_max.max(non_conservation_amount.clone().max());

            for i in 0..3 {
                conservation_check[i].push((momentum[0][n].0, non_conservation_amount[i]));
            }
        }
        y_min = y_min.min(-1e-15);
        y_max = y_max.max(1e-15);
        (conservation_check, (y_min, y_max))
    };

    let colors: [RGBColor; 3] = [BLUE, GREEN, RED];
    let coord_names = vec!["x","y","z"];
            
    let mut momentum_chart = ChartBuilder::on(&areas[0])
        .caption("Momentum", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..max_time, y_min * 1.1..y_max * 1.1)?;

    momentum_chart.configure_mesh()
        .x_desc(format!("Time ({time_unit})"))           
        .y_desc("Momentum (M_sun km / s)")
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    for i in 0..3 {
        let color = colors[i].clone();
        momentum_chart.draw_series(LineSeries::new(momentum[i].clone(), &colors[i]))?.label(format!("p_{}", coord_names[i])).legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color));
    }
        
    momentum_chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    let mut conservation_chart = ChartBuilder::on(&areas[1])
        .caption("Momentum Conservation", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..max_time, conservation_bounds.0 * 1.1..conservation_bounds.1 * 1.1)?;

    dbg!("before");
    conservation_chart.configure_mesh()
        .x_desc(format!("Time ({time_unit})"))           
        .y_desc("(P - P_0) / |P_0|")
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;
    dbg!("after");

    for i in 0..3 {
        conservation_chart.draw_series(LineSeries::new(conservation_check[i].clone(), colors[i]))?;
    }
    
    root.present()?;
    println!("Momentum plot saved as {}", filename);
    Ok(())
}

pub fn plot_angular_momentum(data_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (2*1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let areas = root.split_evenly((1,2));

    let mut angular_momentum: Vec<Vec<(f64, f64)>> = vec![Vec::new(); 3];
    
    let mut file_num = 0; 
    let mut max_time = 0.0;
    let mut y_max = f64::MIN;
    let mut y_min = f64::MAX;
    while Path::new(&format!("{}/snapshot_{:03}.log", data_directory, file_num)).exists() {
        let state = load_checkpoint(data_directory, &file_num)?;

        let angular_momentum_data = calculate_angular_momentum(&state.data);

        for i in 0..3 {
            angular_momentum[i].push((state.time, angular_momentum_data[i]));
        }
        
        y_max = y_max.max(angular_momentum_data.max());
        y_min = y_min.min(angular_momentum_data.min());

        max_time = state.time;
        file_num += 1
    }
    y_max = y_max.max(1e-3);
    y_min = y_min.min(-1e-3);
    

    let mut metric_factor = 0;
    while max_time < (10.0 as f64).powi(metric_factor) {
        metric_factor -= 3
    }
    for n in 0..angular_momentum[0].len() {
        for i in 0..3 {
            angular_momentum[i][n].0 /= (10.0 as f64).powi(metric_factor);
        }
    }
    let time_unit = match metric_factor {
        0 => "Gyr",
        -3 => "Myr",
        -6 => "Kyr",
        -9 => "Year",
        _ => &format!("1e{} Years", 9 - metric_factor) // largest supported unit
    };

    let (conservation_check , conservation_bounds) = {
        let mut conservation_check = vec![Vec::with_capacity(angular_momentum[0].len()), Vec::with_capacity(angular_momentum[0].len()),Vec::with_capacity(angular_momentum[0].len())];
        let mut y_min = f64::MAX;
        let mut y_max = f64::MIN;
        for n in 0..angular_momentum[0].len() {
            let non_conservation_amount = {
                let mut non_conservation_amount = Vec::new();
                for i in 0..3 {
                    non_conservation_amount.push((angular_momentum[i][n].1 - angular_momentum[i][0].1) / angular_momentum[i][0].1.abs().max(1e-10)) //avoid division by 0
                }
                non_conservation_amount
            };

            y_min = y_min.min(non_conservation_amount.clone().min());
            y_max = y_max.max(non_conservation_amount.clone().max());

            for i in 0..3 {
                conservation_check[i].push((angular_momentum[0][n].0, non_conservation_amount[i]));
            }
        }
        y_max = y_max.max(1e-15);
        y_min = y_min.min(-1e-15);
        (conservation_check, (y_min, y_max))
    };

    let colors: [RGBColor; 3] = [BLUE, GREEN, RED];
    let coord_names = vec!["x","y","z"];
            
    let mut angular_momentum_chart = ChartBuilder::on(&areas[0])
        .caption("Angular Momentum", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..max_time, y_min * 1.1..y_max * 1.1)?;

    angular_momentum_chart.configure_mesh()
        .x_desc(format!("Time ({time_unit})"))           
        .y_desc("Angular Momentum (M_sun kpc km / s)")
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    for i in 0..3 {
        let color = colors[i].clone();
        angular_momentum_chart.draw_series(LineSeries::new(angular_momentum[i].clone(), &colors[i]))?.label(format!("L_{}", coord_names[i])).legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], color));
    }
        
    angular_momentum_chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    let mut conservation_chart = ChartBuilder::on(&areas[1])
        .caption("Angular Momentum Conservation", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..max_time, conservation_bounds.0 * 1.1..conservation_bounds.1 * 1.1)?;

    conservation_chart.configure_mesh()
        .x_desc(format!("Time ({time_unit})"))           
        .y_desc("(L - L_0) / |L_0|")
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    for i in 0..3 {
        conservation_chart.draw_series(LineSeries::new(conservation_check[i].clone(), colors[i]))?;
    }
    
    root.present()?;
    println!("Angular Momentum plot saved as {}", filename);
    Ok(())
}

pub fn plot_function(x_points: &Vec<f64>, y_points: &Vec<f64>, filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>> { 
    let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;

    for i in 0..y_points.len() {
        if y_points[i] > y_max {
            y_max = y_points[i]
        }
        if y_points[i] < y_min {
            y_min = y_points[i]
        }
    }

    //println!("y_min = {:.3}, y_max={:.3}", y_min, y_max);

    let x_range = x_points[0]..x_points[x_points.len() - 1];

    let y_range = (y_min+1e-4) * match y_min.signum() {
                1.0 => 0.9,
                -1.0 => 1.1,
                _ => panic!("number has no sign, is probably NaN")
            }..y_max * match y_max.signum() {
                1.0 => 1.1,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN")
            };

    //let y_range = y_range.log_scale();

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(x_range, y_range)?;

    chart.configure_mesh()
        .x_desc(xlabel) // X-axis label
        .y_desc(ylabel) // Y-axis label
        .x_label_formatter(&|x| {
        if x.abs() >= 1000.0 || x.abs() <= 0.1 {
            format!("{:.1e}", x)
        } else {
            format!("{:.1}", x)
        }
        })
        .y_label_formatter(&|y| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;

    let plot_profile: Vec<(f64, f64)> = (0..x_points.len())
        .map(|i| (x_points[i], y_points[i]))
        .collect();

    chart.draw_series(LineSeries::new(plot_profile, &BLUE))?;

    root.present()?;
    println!("Plot saved as {}", filename);
    Ok(())
}

pub fn plot_check_function<T: Fn(f64) -> f64>(x_points: &Vec<f64>, analytic_check: &T, numerical_check: &Vec<f64>, filename: &str, title: &str, xlabel: &str, ylabel: &str) -> Result<(), Box<dyn std::error::Error>> { 
    let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut y_min = f64::MAX;
    let mut y_max = f64::MIN;

    let mut analytic_points: Vec<f64> = vec![0.0; x_points.len()];
    for i in 0..x_points.len() {
        analytic_points[i] = analytic_check(x_points[i]);
    }
    
    for i in 0..analytic_points.len() {
        //finding max
        if analytic_points[i] > y_max {
            y_max = analytic_points[i]
        }
        if numerical_check[i] > y_max {
            y_max = numerical_check[i]
        }

        //finding min
        if analytic_points[i] < y_min {
            y_min = analytic_points[i]
        }
        if numerical_check[i] < y_min {
            y_min = numerical_check[i]
        }
    }

    let x_range = (x_points[0]..x_points[x_points.len() - 1]).log_scale();
    
    let y_range = ((y_min+1e-4) * match y_min.signum() {
            1.0 => 0.9,
            -1.0 => 1.1,
            _ => panic!("number has no sign, is probably NaN")
        }..y_max * match y_max.signum() {
            1.0 => 1.1,
            -1.0 => 0.9,
            _ => panic!("number has no sign, is probably NaN")
        }).log_scale();
    
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(60)
        .build_cartesian_2d(x_range, y_range)?;

    chart.configure_mesh()
        .x_desc(xlabel)           // X-axis label
        .y_desc(ylabel) // Y-axis label
        .x_label_formatter(&|x: &f64| {
            if x.abs() >= 1000.0 || x.abs() <= 0.1 {
                format!("{:.1e}", x)
            } else {
                format!("{:.1}", x)
            }
        })
        .y_label_formatter(&|y: &f64| {
            if y.abs() >= 1000.0 || y.abs() <= 0.1 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;
    
    let plot_profile: Vec<(f64, f64)> = (0..x_points.len())
        .map(|i| (x_points[i], numerical_check[i]))
        .collect();

    chart.draw_series(LineSeries::new(plot_profile, &BLUE))?.label("Numerical Result").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    
    let plot_profile: Vec<(f64, f64)> = (0..x_points.len())
        .map(|i| (x_points[i], analytic_points[i]))
        .collect();

    chart.draw_series(DashedLineSeries::new(plot_profile, 5, 5, (&RED).into()))?.label("Analytic Result").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;
    
    root.present()?;
    println!("Analytic check plot saved as {}", filename);
    Ok(())
}