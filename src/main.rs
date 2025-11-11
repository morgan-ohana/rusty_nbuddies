mod black_hole;
mod forces;
mod time_evol;
mod init_conds;
mod plotting;
mod logging;
mod eddington_inverter;
mod diagnostic;

use std::fs;
use std::path::Path;

use crate::black_hole::BlackHole;
use crate::init_conds::*;
use crate::forces::*;
use crate::time_evol::*;
use crate::plotting::*;
use crate::logging::*;
use crate::diagnostic::*;

const AU: f64 = 4.848136811e-9; // AU in kpc

const BATCHES: usize = 100;
const STEPS: usize = BATCHES*BATCH_SIZE;
const BATCH_SIZE: usize = 10000000;
const DELTA_T: f64 = 0.001*YEAR;
const OUTPUT_DIRECTORY: &str = "output";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //abg_profile_init_conds(&1.0, &3.0, &1.0, &1.0, &1e7, &14.1, &4)?;
    plummer_init_conds(&1.0, &1e8, &4)?;

    //run_simulation()?;

    //plot_black_hole_trajectories(OUTPUT_DIRECTORY, "test_trajectories.png")?;
    //create_comprehensive_plots("test", &data, &DELTA_T)?;
    
    //check_energy_conservation(OUTPUT_DIRECTORY)?;
    //check_angular_momentum_conservation(OUTPUT_DIRECTORY)?;

    Ok(())
}

fn run_simulation() -> Result<(), Box<dyn std::error::Error>> {
    clean_output()?;
    set_up_output_directory()?;
    println!("Timestep is {}", DELTA_T);

    //let mut black_holes: [BlackHole; 2] = binary_init_conds(10.0*AU, 3.0e-9, 1.0, 2.0);
    let mut black_holes: Vec<BlackHole> = binary_circular_init_conds(1000.0*AU, 10000.0, 10000.0);
    let mut previous_black_holes: Vec<BlackHole> = black_holes.clone();

    let mut data: Vec<Vec<BlackHole>> = vec![black_holes.clone(); BATCH_SIZE];

    print_all_info(&black_holes);

    for n in 0..STEPS {
        
        data[n % BATCH_SIZE] = black_holes.clone();

        recalculate_acceleration_due_to_gravity(&mut black_holes);

        update_velocities(&mut black_holes, &previous_black_holes, &DELTA_T);

        update_positions(&mut black_holes, &previous_black_holes, &DELTA_T);
        
        previous_black_holes = black_holes.clone();

        if n % BATCH_SIZE == BATCH_SIZE - 1 {
            let sim = SimulationState {
                time: n as f64 * DELTA_T,
                data: data.clone(),
                step_count: n
            };
            save_checkpoint(&sim, &OUTPUT_DIRECTORY, &BATCH_SIZE)?;
        }
    }
    print_all_info(&black_holes);

    Ok(())
}

fn clean_output() -> std::io::Result<()> {
    if Path::new(OUTPUT_DIRECTORY).exists() {
        fs::remove_dir_all(OUTPUT_DIRECTORY)?;
    }
    Ok(())
}

fn set_up_output_directory() -> std::io::Result<()> {
    if !Path::new(OUTPUT_DIRECTORY).exists() {
        fs::create_dir_all(OUTPUT_DIRECTORY)?;
    }
    Ok(())
}

fn print_all_info(black_holes: &Vec<BlackHole>) { 
    println!("x pos {}", black_holes[0].position[0]);
    println!("x vel {}", black_holes[0].velocity[0]);
    println!("x accel {}", black_holes[0].acceleration[0]);
    println!("x pos {}", black_holes[1].position[0]);
    println!("x vel {}", black_holes[1].velocity[0]);
    println!("x accel {}", black_holes[1].acceleration[0]);
}