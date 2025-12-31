use std::thread;
use std::thread::JoinHandle;

use crate::particle::Particle;
use crate::forces::{recalculate_dynamics_due_to_gravity, ForceCalculationMethod};
use crate::logging::{save_checkpoint, SimulationState};
use crate::time_evol::{compute_timestep, update_positions, update_velocities};

pub const GYR: f64 = 3.1536e16; // Myr in seconds

pub fn run_simulation(init_conds: Vec<Particle>, max_time: f64, batch_duration: f64, method: ForceCalculationMethod, eta: f64, output_directory: String) {
    
    let mut particles: Vec<Particle> = init_conds;
    
    //initial dynamics calculation
    recalculate_dynamics_due_to_gravity(&mut particles, &method);
    recalculate_dynamics_due_to_gravity(&mut particles, &method); //two calls to "warm up" higher order derivatives

    let mut running_time: f64 = 0.0;
    let mut step: usize = 0;
    let mut batch_num: usize = 0;
    let mut timestep;
    let mut logging_handles: Vec<JoinHandle<()>> = Vec::new();
    let mut should_log = true;

    while running_time < max_time {

        if should_log {
            let sim = SimulationState {
                time: running_time / GYR,
                data: particles.clone(),
                step_count: step
            };
            
            let output_directory_copy = output_directory.clone();
            logging_handles.push(thread::spawn(move || {
                save_checkpoint(&sim, &output_directory_copy, &batch_num).expect(&format!("Failed to save checkpoint number {batch_num}"));
            }));

            should_log = false;
            batch_num += 1;
        }
        
        // dynamical timestep calculation
        timestep = f64::MAX;
        for m in 0..particles.len() {
            timestep = timestep.min(compute_timestep(&particles[m], &eta));
        }
        
        if running_time + timestep > batch_duration * (batch_num as f64) {
            timestep = batch_duration * (batch_num as f64) - running_time;
            should_log = true;
        }
        running_time += timestep;

        // initial kick v_i -> v_{i+0.5}
        update_velocities(&mut particles, &(0.5 * timestep));
        
        // drift x_i -> x_{i+1}
        update_positions(&mut particles, &timestep);
        
        // update dynamics a_{i+1} = A(x_{i+1})
        recalculate_dynamics_due_to_gravity(&mut particles, &method);
        
        // final kick v_{i+0.5} -> v_{i+1}
        update_velocities(&mut particles, &(0.5 * timestep));

        // increment step counter
        step += 1;
    }

    // log final step
    let sim = SimulationState {
        time: running_time / GYR,
        data: particles.clone(),
        step_count: step
    };
    
    save_checkpoint(&sim, &output_directory, &batch_num).expect(&format!("Failed to save checkpoint number {batch_num}"));
            
    for handle in logging_handles {
        handle.join().expect("Logging thread panicked");
    }
}