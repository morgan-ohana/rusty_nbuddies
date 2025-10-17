use std::path::Path;

use crate::black_hole::BlackHole;
use crate::forces::{calculate_energy, calculate_angular_momentum};
use crate::logging::load_checkpoint;
use crate::plotting::YEAR;

pub fn check_energy_conservation(output_directory: &str) -> anyhow::Result<()> {
    let initial_state = load_checkpoint(output_directory, &0)?;
    let initial_energies = calculate_energy(&initial_state.data[0]);
    let init_energy = initial_energies.0 + initial_energies.1;

    let mut last_file_num = 0;
    while Path::new(&format!("{}/restart_{:03}.log", output_directory, last_file_num+1)).exists() {
        last_file_num += 1;
    }

    let final_state = load_checkpoint(output_directory, &last_file_num)?;
    let final_energies = calculate_energy(&final_state.data.last().expect("Somehow the data of the last file is empty"));
    let final_energy = final_energies.0 + final_energies.1;

    let energy_percent_error = (final_energy - init_energy) / f64::abs(init_energy);
    println!("percent energy error per year: {}", energy_percent_error / (final_state.time / YEAR));
    println!("total percent energy error over {:03} years: {}", (final_state.time / YEAR), energy_percent_error);

    Ok(())
}

pub fn check_angular_momentum_conservation(output_directory: &str) -> anyhow::Result<()> {
    let initial_state = load_checkpoint(output_directory, &0)?;
    let initial_angular_momentum = calculate_angular_momentum(&initial_state.data[0]);

    let mut last_file_num = 0;
    while Path::new(&format!("{}/restart_{:03}.log", output_directory, last_file_num+1)).exists() {
        last_file_num += 1;
    }

    let final_state = load_checkpoint(output_directory, &last_file_num)?;
    let final_angular_momentum = calculate_angular_momentum(&final_state.data.last().expect("Somehow the data of the last file is empty"));
    
    let mut angular_momentum_percent_error: [f64; 3] = [0.0, 0.0, 0.0];
    for i in 0..3 {
        angular_momentum_percent_error[i] = (final_angular_momentum[i] - initial_angular_momentum[i])/ f64::abs(initial_angular_momentum[i]);
    }

    let ang_mom_names = ["L_x", "L_y", "L_z"];

    for i in 0..3 {
        println!("percent {} error per year: {}", ang_mom_names[i], angular_momentum_percent_error[i] / (final_state.time / YEAR));
        println!("total percent {} error over {:03} years: {}", ang_mom_names[i], (final_state.time / YEAR), angular_momentum_percent_error[i]);
    }
    
    Ok(())
}