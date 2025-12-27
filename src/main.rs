mod particle;
mod gravitree;
mod forces;
mod time_evol;
mod init_conds;
mod plotting;
mod logging;
mod eddington_inverter;
mod diagnostic;
mod vectors;

use std::fs;
use std::path::Path;
use std::thread;
use std::thread::JoinHandle;
use std::time::Instant;

use crate::particle::Particle;
use crate::gravitree::{AccuracyCriterion, Node, };
use crate::init_conds::*;
use crate::forces::*;
use crate::time_evol::*;
use crate::plotting::*;
use crate::logging::*;
use crate::diagnostic::*;

const AU: f64 = 4.848136811e-9; // AU in kpc

const ETA: f64 = 0.01; //timestep accuracy parameter
const BATCHES: usize = 100;
const STEPS: usize = BATCHES*BATCH_SIZE;
const BATCH_SIZE: usize = 100;
const OUTPUT_DIRECTORY: &str = "output";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //nfw_init_conds(1.0, 1e7, 14.1, 1000000);
    //plummer_init_conds(1.0, 1e8, 1000000);
    //abg_profile_init_conds(&15.0, &4.5, &2.0, &3.14, &4e9, Some(&5.3), &2);

    //let start = Instant::now();
    //run_simulation()?;
    //let duration = start.elapsed();
    //println!("Time elapsed in simulation() is: {:?}", duration);

    plot_trajectories(OUTPUT_DIRECTORY, "test_trajectories.png")?;
    //create_comprehensive_plots("test", OUTPUT_DIRECTORY)?;
    
    //check_energy_conservation(OUTPUT_DIRECTORY)?;
    //check_angular_momentum_conservation(OUTPUT_DIRECTORY)?;

    Ok(())
}

fn run_simulation() -> Result<(), Box<dyn std::error::Error>> {
    clean_output()?;
    set_up_output_directory()?;

    //let method = ForceCalculationMethod::Direct;
    let method = ForceCalculationMethod::Tree(AccuracyCriterion::Dynamical(1e-3));

    //let mut particles: Vec<Particle> = binary_init_conds(10.0*AU, 3.0e-9, 1.0, 2.0);
    //let mut particles: Vec<Particle> = binary_circular_init_conds(1.0, 1e7, 1e7);
    //let mut particles: Vec<Particle> = binary_circular_init_conds(10.0, 1e8, 1e8);
    let mut particles: Vec<Particle> = plummer_init_conds(1.0, 1e8, 10000);
    //let mut particles: Vec<Particle> = nfw_init_conds(1.0, 1e7, 14.1, 100);

    //print_all_info(&particles);
    //initial dynamics calculation
    recalculate_dynamics_due_to_gravity(&mut particles, &method);
    recalculate_dynamics_due_to_gravity(&mut particles, &method); //two calls to "warm up" higher order derivatives

    let mut running_time: f64 = 0.0;
    let mut timestep;
    let mut logging_handles: Vec<JoinHandle<()>> = Vec::new();

    for n in 0..STEPS {

        if n % BATCH_SIZE == 0 {
            let sim = SimulationState {
                time: running_time,
                data: particles.clone(),
                step_count: n
            };
            logging_handles.push(thread::spawn(move || {
                save_checkpoint(&sim, &OUTPUT_DIRECTORY, &BATCH_SIZE).expect("Failed to save checkpoint");
            }));
        }

        timestep = f64::MAX;
        for m in 0..particles.len() {
            timestep = timestep.min(compute_timestep(&particles[m], &ETA))
        }
        running_time += timestep;

        //initial kick v_i -> v_{i+0.5}
        update_velocities(&mut particles, &(0.5 * timestep));
        
        //drift x_i -> x_{i+1}
        update_positions(&mut particles, &timestep);
        
        //update dynamics a_{i+1} = A(x_{i+1})
        recalculate_dynamics_due_to_gravity(&mut particles, &method);
        
        //final kick v_{i+0.5} -> v_{i+1}
        update_velocities(&mut particles, &(0.5 * timestep));

    }

    for handle in logging_handles {
        handle.join().expect("Logging thread panicked");
    }
    //print_all_info(&particles);

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

fn print_all_info(particles: &Vec<Particle>) { 
    println!("--- Particle Info ---");
    for particle in particles {
        println!("  x pos {}", particle.position[0]);
        println!("  x vel {}", particle.velocity[0]);
        println!("  x accel {}", particle.acceleration[0]);
    }
}