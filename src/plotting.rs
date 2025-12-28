use plotters::prelude::*;
use std::path::Path;
use std::f64::consts::PI;
use crate::particle;
use crate::particle::Particle;
use crate::forces::calculate_energy;
use crate::forces::GG;
use crate::logging::load_checkpoint;

pub const YEAR: f64 = 31556952.0; // 1 year in sec

pub fn plot_trajectories(output_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let scale = 4;
    let root = BitMapBackend::new(filename, (1024*scale, 1024*scale)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut trajectories: Vec<Vec<(f64, f64, f64)>> = Vec::new();
    let mut bounds: (f64, f64) = (0.0, 0.0);
    
    let mut i: usize = 0;
    let mut last_data: Vec<Particle> = Vec::new();

    while Path::new(&format!("{}/snapshot_{:03}.log", output_directory, i)).exists() {
        println!("{}", i);
        let data = load_checkpoint(output_directory, &i)?.data;

        if trajectories.is_empty() {
            trajectories = vec![Vec::new(); data.len()];
        }

        // Add current positions to each particle's trajectory
        for (particle_idx, particle) in data.iter().enumerate() {
            trajectories[particle_idx].push((
                particle.position[0],
                particle.position[1],
                particle.position[2],
            ));

            bounds.0 = bounds.0.min(particle.position[0]).min(particle.position[1]);
            bounds.1 = bounds.1.max(particle.position[0]).max(particle.position[1]);
        }
        
        last_data = data;
        i += 1;
    }

    let mut chart = ChartBuilder::on(&root)
        .caption("Black Hole Binary System Trajectories", ("sans-serif", 40))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_3d(bounds.0..bounds.1, bounds.0..bounds.1, bounds.0..bounds.1)?;

    chart.configure_axes().draw()?;

    let colors = vec![&RED, &BLUE, &GREEN, &YELLOW, &CYAN, &MAGENTA];

    //let particle_count = trajectories.len();

    // Plot trajectories for each particle
    for (particle_idx, traj) in trajectories.iter().enumerate() {
        if traj.is_empty() {
            continue;
        }
        
        let color = colors[particle_idx % colors.len()];
        
        chart.draw_series(LineSeries::new(
            traj.clone(),
            color.stroke_width(1),
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
            2,
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

pub fn create_comprehensive_plots(name: &str, output_directory: &str) -> Result<(), Box<dyn std::error::Error>> {
    // 1. Trajectory plot
    //plot_black_hole_trajectories(output_directory, &(name.to_owned() + "_trajectories.png"))?;
    
    // 2. Distance between black holes over time
    //plot_separation_vs_time(output_directory, &(name.to_owned() + "_separation_vs_time.png"))?;

    // 3. Angular Momentum conservation check
    //plot_angular_momentum_vs_time(output_directory, &(name.to_owned() + "_angular_momentum_vs_time.png"))?;
    
    // 4. Energy conservation check
    //plot_energy_vs_time(output_directory, &(name.to_owned() + "_energy_vs_time.png"))?;
    
    Ok(())
}

// pub fn plot_separation_vs_time(output_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
//     let state = load_checkpoint(output_directory, &0)?;
//     let data = state.data;
    
//     if data.len() != 2 {
//         panic!("Separation plot only implemented for 2-body systems");
//     }

//     let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
//     root.fill(&WHITE)?;

//     let separations: Vec<f64> = data.iter()
//         .map(|frame| {
//             let dx = frame[0].position[0] - frame[1].position[0];
//             let dy = frame[0].position[1] - frame[1].position[1];
//             (dx * dx + dy * dy).sqrt()
//         })
//         .collect();
    
//     let max_time = (data.len() as f64) * state.time / YEAR;
//     let max_sep = separations.iter().cloned().fold(f64::MIN, f64::max);
//     let min_sep = separations.iter().cloned().fold(f64::MAX, f64::min);

//     let mut chart = ChartBuilder::on(&root)
//         .caption("Separation Between Black Holes", ("sans-serif", 40))
//         .margin(10)
//         .x_label_area_size(30)
//         .y_label_area_size(60)
//         .build_cartesian_2d(0.0..max_time, min_sep * 0.9..max_sep * 1.1)?;

//     chart.configure_mesh()
//         .x_desc("Time (years)")           // X-axis label
//         .y_desc("Separation Distance (kpc) ") // Y-axis label
//         .x_label_formatter(&|x| {
//         if x.abs() >= 1000.0 {
//             format!("{:.1e}", x)
//         } else {
//             format!("{:.1}", x)
//         }
//         })
//         .y_label_formatter(&|y| {
//             if y.abs() >= 1000.0 {
//                 format!("{:.1e}", y)
//             } else {
//                 format!("{:.1}", y)
//             }
//         })
//         .draw()?;

//     let mut time_points: Vec<(f64, f64)> = (0..data.len())
//         .map(|i| (i as f64, separations[i]))
//         .collect();

//     for i in 0..data.len() {
//         time_points[i].0 *= state.time / YEAR;
//     }

//     chart.draw_series(LineSeries::new(time_points, &BLUE))?;

//     root.present()?;
//     println!("Separation plot saved as {}", filename);
//     Ok(())
// }

// pub fn plot_angular_momentum_vs_time(output_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
//     let state = load_checkpoint(output_directory, &0)?;
//     let data = state.data;
    
//     let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
//     root.fill(&WHITE)?;

//     let angular_momentum: Vec<f64> = data.iter()
//         .map(|frame| {
//             let l_1 = frame[0].position[0] * frame[0].mass * frame[0].velocity[1] - frame[0].position[1] * frame[0].mass * frame[0].velocity[0];
//             let l_2 = frame[1].position[0] * frame[1].mass * frame[1].velocity[1] - frame[1].position[1] * frame[1].mass * frame[1].velocity[0];
//             l_1 + l_2
//         })
//         .collect();
    
//     let max_time = (data.len() as f64) * state.time / YEAR;
//     let max_ang_mom = angular_momentum.iter().cloned().fold(f64::MIN, f64::max);
//     let min_ang_mom = angular_momentum.iter().cloned().fold(f64::MAX, f64::min);

//     let mut chart = ChartBuilder::on(&root)
//         .caption("Total Angular Momentum Black Holes", ("sans-serif", 40))
//         .margin(10)
//         .x_label_area_size(30)
//         .y_label_area_size(60)
//         .build_cartesian_2d(0.0..max_time, min_ang_mom * 0.9..max_ang_mom * 1.1)?;

//     chart.configure_mesh()
//         .x_desc("Time (years)")           // X-axis label
//         .y_desc("Angular Momentum (M_sun km kpc / s) ") // Y-axis label
//         .x_label_formatter(&|x| {
//         if x.abs() >= 1000.0 {
//             format!("{:.1e}", x)
//         } else {
//             format!("{:.1}", x)
//         }
//         })
//         .y_label_formatter(&|y| {
//             if y.abs() >= 1000.0 {
//                 format!("{:.1e}", y)
//             } else {
//                 format!("{:.1}", y)
//             }
//         })
//         .draw()?;

//     let mut time_points: Vec<(f64, f64)> = (0..data.len())
//         .map(|i| (i as f64, angular_momentum[i]))
//         .collect();

//     for i in 0..data.len() {
//         time_points[i].0 *= state.time / YEAR;
//     }

//     chart.draw_series(LineSeries::new(time_points, &BLUE))?;

//     root.present()?;
//     println!("Angular Momentum plot saved as {}", filename);
//     Ok(())
// }

// pub fn plot_energy_vs_time(output_directory: &str, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
//     let state = load_checkpoint(output_directory, &0)?;
//     let data = state.data;
    
//     let root = BitMapBackend::new(filename, (1024, 768)).into_drawing_area();
//     root.fill(&WHITE)?;

//     let kinetic: Vec<f64> = data.iter()
//         .map(|frame| {
//             calculate_energy(frame).0
//         })
//         .collect();

//     let potential: Vec<f64> = data.iter()
//         .map(|frame| {
//             calculate_energy(frame).1
//         })
//         .collect();

//     let total_energy: Vec<f64> = data.iter()
//         .map(|frame| {
//             calculate_energy(frame).0 + calculate_energy(frame).1
//         })
//         .collect();
    
//     let max_time = (data.len() as f64) * state.time / YEAR;
//     let ymax = kinetic.iter().cloned().fold(f64::MIN, f64::max);
//     let ymin = potential.iter().cloned().fold(f64::MAX, f64::min);

//     let mut chart = ChartBuilder::on(&root)
//         .caption("Total Energy of Black Holes", ("sans-serif", 40))
//         .margin(10)
//         .x_label_area_size(30)
//         .y_label_area_size(60)
//         .build_cartesian_2d(0.0..max_time, ymin * 0.9..ymax * 1.1)?;

//     chart.configure_mesh()
//         .x_desc("Time (years)")           
//         .y_desc("Energy (M_sun km^2 / s^2)")
//         .x_label_formatter(&|x| {
//         if x.abs() >= 1000.0 {
//             format!("{:.1e}", x)
//         } else {
//             format!("{:.1}", x)
//         }
//         })
//         .y_label_formatter(&|y| {
//             if y.abs() >= 1000.0 {
//                 format!("{:.1e}", y)
//             } else {
//                 format!("{:.1}", y)
//             }
//         })
//         .draw()?;

//     let mut kinetic_points: Vec<(f64, f64)> = (0..data.len())
//         .map(|i| (i as f64, kinetic[i]))
//         .collect();

//     let mut potential_points: Vec<(f64, f64)> = (0..data.len())
//         .map(|i| (i as f64, potential[i]))
//         .collect();

//     let mut energy_points: Vec<(f64, f64)> = (0..data.len())
//         .map(|i| (i as f64, total_energy[i]))
//         .collect();

//     println!("Delta E/Year = {}", (energy_points.last().expect("Energy series empty in energy plot").1 - energy_points[0].1)/(energy_points.last().expect("Energy series empty in energy plot").0 - energy_points[0].0));

//     for i in 0..data.len() {
//         kinetic_points[i].0 *= state.time / YEAR;
//         potential_points[i].0 *= state.time / YEAR;
//         energy_points[i].0 *= state.time / YEAR;
//     }

//     chart.draw_series(LineSeries::new(kinetic_points, &BLUE))?.label("Kinetic").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
//     chart.draw_series(LineSeries::new(potential_points, &GREEN))?.label("Potential").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
//     chart.draw_series(LineSeries::new(energy_points, &BLACK))?.label("Total").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

//     chart.configure_series_labels()
//         .background_style(&WHITE.mix(0.8))
//         .border_style(&BLACK)
//         .draw()?;

//     root.present()?;
//     println!("Energy plot saved as {}", filename);
//     Ok(())
// }

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

    let x_range = (x_points[0]..x_points[x_points.len() - 1]);

    let y_range = ((y_min+1e-4) * match y_min.signum() {
                1.0 => 0.9,
                -1.0 => 1.1,
                _ => panic!("number has no sign, is probably NaN")
            }..y_max * match y_max.signum() {
                1.0 => 1.1,
                -1.0 => 0.9,
                _ => panic!("number has no sign, is probably NaN")
            });//.log_scale();

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
            if x.abs() >= 1000.0 {
                format!("{:.1e}", x)
            } else {
                format!("{:.1}", x)
            }
        })
        .y_label_formatter(&|y: &f64| {
            if y.abs() >= 1000.0 {
                format!("{:.1e}", y)
            } else {
                format!("{:.1}", y)
            }
        })
        .draw()?;
    
    let mut plot_profile: Vec<(f64, f64)> = (0..x_points.len())
        .map(|i| (x_points[i], numerical_check[i]))
        .collect();

    chart.draw_series(LineSeries::new(plot_profile, &BLUE))?.label("Numerical Result").legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    
    let mut plot_profile: Vec<(f64, f64)> = (0..x_points.len())
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