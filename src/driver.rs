use std::thread;
use crossbeam::channel::{bounded};
use std::sync::Arc;

use crate::particle::Particle;
use crate::gravitree::build_gravitree;
use crate::forces::{ForceCalculationMethod, recalculate_dynamics_due_to_gravity_directly, recalculate_dynamics_due_to_gravity_with_tree};
use crate::logging::{save_checkpoint, SimulationState};
use crate::time_evol::{compute_timestep, update_positions, update_velocities};

pub const GYR: f64 = 3.1536e16; // Myr in seconds

pub fn run_simulation(init_conds: Vec<Particle>, max_time: f64, batch_duration: f64, method: ForceCalculationMethod, eta: f64, output_directory: String) {
    
    let mut particles: Vec<Particle> = init_conds;
    
    //initial dynamics calculation
    //two calls to "warm up" higher order derivatives
    match method {
        ForceCalculationMethod::Direct => {
            recalculate_dynamics_due_to_gravity_directly(&mut particles);
            recalculate_dynamics_due_to_gravity_directly(&mut particles);
        },
        ForceCalculationMethod::Tree(ref criterion) => {
            let root = Some(Box::new(build_gravitree(particles.clone())));
            recalculate_dynamics_due_to_gravity_with_tree(&mut particles, &criterion, &root);
            let root = Some(Box::new(build_gravitree(particles.clone())));
            recalculate_dynamics_due_to_gravity_with_tree(&mut particles, &criterion, &root);
        }
    }
    
    let mut running_time: f64 = 0.0;
    let mut step: usize = 0;
    let mut batch_num: usize = 0;
    let mut timestep;
    let mut should_log = true;

    // Create a channel for sending work to logger threads
    // The buffer size controls how many checkpoints can queue up
    let (tx, rx) = bounded(10); // Buffer 10 checkpoints

    // Define thread pool size
    let num_loggers = 1;
    let output_dir_arc = Arc::new(output_directory.clone());
    let mut logger_handles = Vec::new();
    
    // Spawn logger threads
    for logger_id in 0..num_loggers {
        let rx = rx.clone(); // Each logger gets its own receiver
        let output_dir = output_dir_arc.clone();
        
        logger_handles.push(thread::spawn(move || {
            // Each logger runs this loop
            while let Ok((sim, batch_num)) = rx.recv() {
                println!("Logger {} saving checkpoint {}", logger_id, batch_num);
                save_checkpoint(&sim, &output_dir, &batch_num)
                    .unwrap_or_else(|e| eprintln!("Failed to save checkpoint {}: {}", batch_num, e));
            }
            println!("Logger {} shutting down", logger_id);
        }));
    }
    
    // Drop the original receiver (we're done cloning it for loggers)
    drop(rx); // Important: This ensures loggers exit when we drop tx later

    while running_time < max_time {

        if should_log {
            let sim = SimulationState {
                time: running_time / GYR,
                data: particles.clone(),
                step_count: step
            };
            
            // Send to logger thread instead of spawning new thread
            tx.send((sim, batch_num)).expect("Failed to send checkpoint to logger");

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
        match method {
            ForceCalculationMethod::Direct => recalculate_dynamics_due_to_gravity_directly(&mut particles),
            ForceCalculationMethod::Tree(ref criterion) => {
                let root = Some(Box::new(build_gravitree(particles.clone())));
                recalculate_dynamics_due_to_gravity_with_tree(&mut particles, &criterion, &root);
            }
        }
        
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
            
    // Signal logger to stop and wait for it
    drop(tx); // This will cause rx.recv() to return Err, breaking the logger loop
    for handle in logger_handles {
        handle.join().expect("logger thread panicked");
    }
}