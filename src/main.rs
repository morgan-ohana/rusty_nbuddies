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

use core::time;
use std::fs;
use std::path::Path;
use std::process::Output;
use std::thread;
use std::thread::JoinHandle;
use std::time::Instant;
use std::env;
use config::Config;
use statrs::distribution::Geometric;
//use std::collections::HashMap;

use crate::particle::Particle;
use crate::gravitree::{AccuracyCriterion, Node, };
use crate::init_conds::*;
use crate::forces::*;
use crate::time_evol::*;
use crate::plotting::*;
use crate::logging::*;
use crate::diagnostic::*;

const AU: f64 = 4.848136811e-9; // AU in kpc
const GYR: f64 = 3.1536e16; // Myr in seconds

// Default config values
const DEFAULT_ETA: f64 = 0.01; // timestep accuracy parameter
const DEFAULT_ALPHA: f64 = 1e-3; // dynamical tree accuracy parameter
const DEFAULT_THETA: f64 = 0.3; // geometric tree accuracy parameter 
const DEFAULT_BATCH_DURATION: f64 = 0.1 * GYR; // time per batch
const DEFAULT_MAX_TIME: f64 = 30.0 * GYR;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let config_path = &args[1];

    let settings = Config::builder()
        // Add in `./Settings.toml`
        .add_source(config::File::with_name(config_path))
        // Add in settings from the environment (with a prefix of APP)
        // Eg.. `APP_DEBUG=1 ./target/app` would set the `debug` key
        .add_source(config::Environment::with_prefix("APP"))
        .build()
        .unwrap();

    // Prepare Initial Conditions
    let init_conds_table = settings.get_table("Initial-Conditions").expect("Initial conditions table must be present in config file, please consult the documentation for formatting guidelines.");
    let init_conds_file: String = match init_conds_table.get("type")
                                        .expect("No type found for initial conditions, ensure you include type = {valid type string} in your initial conditions table")
                                        .clone().into_string().expect("Initial conditions type must be a string, consult documentation for valid types.").as_str() {
        "File" => {
            init_conds_table.get("location")
                .expect("Unable to read file location, ensure you include location = {some string} in your initial conditions table")
                .clone().into_string().expect("Location of initial conditions file must be a string!")
        },
        "Plummer" => {
            let r_s = init_conds_table.get("r_s")
                    .expect("Unable to read scale radius, ensure you include r_s = {some float} in your initial conditions table")
                    .clone().into_float().expect("Scale radius must be a float!");
            let total_mass = init_conds_table.get("total_mass")
                    .expect("Unable to read total mass, ensure you include total_mass = {some float} in your initial conditions table")
                    .clone().into_float()
                    .expect("Total mass must be a float!");
            let particle_num  = init_conds_table.get("particle_num")
                    .expect("Unable to read particle number, ensure you include particle_num = {some positive int} in your initial conditions table")
                    .clone().into_uint().expect("Particle number must be a positive integer!");
            let output_path = init_conds_table.get("output_directory")
                    .expect("Ubable to read output path, ensure you include output_directory = {some string} in your inital conditions table")
                    .clone().into_string().expect("Output path must be a string!");
            let particles = plummer_init_conds(r_s, total_mass, particle_num as usize, output_path.clone());
            let file_name = output_path.clone() + "/Plummer.log";
            save_init_conds(file_name.clone(), particles).unwrap();
            file_name
        },
        "NFW" => {
            let r_s = init_conds_table.get("r_s")
                    .expect("Unable to read scale radius, ensure you include r_s = {some float} in your initial conditions table")
                    .clone().into_float().expect("Scale radius must be a float!");
            let rho_s = init_conds_table.get("rho_s")
                    .expect("Unable to read scale density, ensure you include rho_s = {some float} in your initial conditions table")
                    .clone().into_float().expect("Scale density must be a float!");
            let r_cutoff = init_conds_table.get("r_cutoff")
                    .expect("Unable to read cutoff radius, ensure you include r_cutoff = {some float} in your initial conditions table")
                    .clone().into_float().expect("Cutoff radius must be a float!");
            let particle_num  = init_conds_table.get("particle_num")
                    .expect("Unable to read particle number, ensure you include particle_num = {some positive int} in your initial conditions table")
                    .clone().into_uint().expect("Particle number must be a positive integer!");
            let output_path = init_conds_table.get("output_directory")
                    .expect("Ubable to read output path, ensure you include output_directory = {some string} in your inital conditions table")
                    .clone().into_string().expect("Output path must be a string!");
            let particles = nfw_init_conds(r_s, rho_s, r_cutoff, particle_num as usize, output_path.clone());
            let file_name = output_path + "/NFW.log";
            save_init_conds(file_name.clone(), particles).unwrap();
            file_name
        },
        "Alpha-Beta-Gamma" => {
            let r_s = init_conds_table.get("r_s")
                    .expect("Unable to read scale radius, ensure you include r_s = {some float} in your initial conditions table")
                    .clone().into_float().expect("Scale radius must be a float!");
            let rho_s = init_conds_table.get("rho_s")
                    .expect("Unable to read scale density, ensure you include rho_s = {some float} in your initial conditions table")
                    .clone().into_float().expect("Scale density must be a float!");
            let alpha = init_conds_table.get("alpha")
                    .expect("Unable to read alpha, ensure you include alpha = {some float} in your initial conditions table")
                    .clone().into_float().expect("Alpha must be a float!");
            let beta = init_conds_table.get("beta")
                    .expect("Unable to read beta, ensure you include beta = {some float} in your initial conditions table")
                    .clone().into_float().expect("Beta must be a float!");
            let gamma = init_conds_table.get("gamma")
                    .expect("Unable to read gamma, ensure you include gamma = {some float} in your initial conditions table")
                    .clone().into_float().expect("Gamma must be a float!");
            let r_cutoff_option = match init_conds_table.get("r_cutoff") {
                Some(r_cutoff) => {
                    Some(r_cutoff.clone().into_float().expect("Cutoff radius must be a float!"))
                },
                None => None
            };
            let particle_num  = init_conds_table.get("particle_num")
                    .expect("Unable to read particle number, ensure you include particle_num = {some positive int} in your initial conditions table")
                    .clone().into_uint().expect("Particle number must be a positive integer!");
            let output_path = init_conds_table.get("output_directory")
                    .expect("Ubable to read output path, ensure you include output_directory = {some string} in your inital conditions table")
                    .clone().into_string().expect("Output path must be a string!");
            let particles = abg_profile_init_conds(alpha, beta, gamma, r_s, rho_s, r_cutoff_option, particle_num as usize, output_path.clone());
            let file_name = output_path + "/Alpha-Beta-Gamma.log";
            save_init_conds(file_name.clone(), particles).unwrap();
            file_name
        },
        _ => {
            panic!("Unknown initial condition type, please consult the documentation for valid initial condition specification.")        
        }
    };
    
    println!("Loading initial conditions from {init_conds_file}");

    let init_conds_data: Vec<Particle> = load_file(init_conds_file).expect("Failed to load initial conditions file, are you sure the path is correct?").data;

    if let Ok(simulation_table) = settings.get_table( "Simulation") {
        let method: ForceCalculationMethod = match simulation_table.get("calculation-method") {
            Some(method_name) => {
                match method_name.clone().into_string().expect("calculation-method must be a string!").as_str() {
                    "Direct" => ForceCalculationMethod::Direct,
                    "Tree" => {
                        match simulation_table.get("accuracy-criterion") {
                            Some(accuracy_criterion) => {
                                match accuracy_criterion.clone().into_string().expect("Accuracy criterion must be a string!").as_str() {
                                    "Geometric" => {
                                        match simulation_table.get("tree-accuracy-parameter") {
                                            Some(accuracy_parameter) => ForceCalculationMethod::Tree(AccuracyCriterion::Geometric(accuracy_parameter.clone().into_float().expect("Geometric accuracy parameter theta must be a float!"))),
                                            // Default parameter for geometric tree
                                            None => ForceCalculationMethod::Tree(AccuracyCriterion::Geometric(DEFAULT_THETA))
                                        }
                                    },
                                    "Dynamical" => {
                                        match simulation_table.get("tree-accuracy-parameter") {
                                            Some(accuracy_parameter) => ForceCalculationMethod::Tree(AccuracyCriterion::Dynamical(accuracy_parameter.clone().into_float().expect("Dynamical accuracy parameter theta must be a float!"))),
                                            // Default parameter for dynamical tree
                                            None => ForceCalculationMethod::Tree(AccuracyCriterion::Dynamical(DEFAULT_ALPHA))
                                        }
                                    },
                                    _ => panic!("Unrecognized accuracy criterion, please consult documentation for valid accuracy criteria.")
                                }
                            },
                            // Default behavior with Tree specified
                            None => ForceCalculationMethod::Tree(AccuracyCriterion::Dynamical(1e-3))
                        }
                    },
                    _ => panic!("Unrecognized calculation method, please consult documentation for valid calculation methods.")
                }
            },
            // Default Behavior
            None => ForceCalculationMethod::Tree(AccuracyCriterion::Dynamical(1e-3))
        };

        let eta = match simulation_table.get("timestep-accuracy-parameter") {
            Some(eta) => eta.clone().into_float().expect("timestep accuracy parameter eta must be a float!"),
            None => DEFAULT_ETA
        };

        let batch_duration = match simulation_table.get("batch-duration") {
            Some(duration) => duration.clone().into_float().expect("batch duration must be a float!"),
            None => DEFAULT_BATCH_DURATION
        };

        let max_time = match simulation_table.get("max-time") {
            Some(max_time) => max_time.clone().into_float().expect("Max time must be a float!") * GYR,
            None => DEFAULT_MAX_TIME
        };

        let output_directory = simulation_table.get("output_directory")
            .expect("Ubable to read output path, ensure you include output_directory = {some string} in your simulation table")
            .clone().into_string().expect("Output path must be a string!");

        let mut startup_message = String::from("Begining simulation:");
        startup_message += &(String::from(format!("\n\t output_directory {output_directory}")));
        startup_message += &(String::from(format!("\n\t Calculation Method: ") + &method.name()));
        startup_message += &(String::from(format!("\n\t Timestep Accuracy Parameter {eta}")));
        startup_message += &(String::from(format!("\n\t Batch Duration {}", batch_duration/GYR)));
        startup_message += &(String::from(format!("\n\t Max Time {}", max_time/GYR)));

        println!("{startup_message}");

        run_simulation(init_conds_data, max_time, batch_duration, method, eta, output_directory)?;
    }

    plot_trajectories("/home/mohana/r_nbud_usr_test/output", "test_trajectories.png")?;
    //create_comprehensive_plots("test", OUTPUT_DIRECTORY)?;
    
    //check_energy_conservation(OUTPUT_DIRECTORY)?;
    //check_angular_momentum_conservation(OUTPUT_DIRECTORY)?;

    Ok(())
}

fn run_simulation(init_conds: Vec<Particle>, max_time: f64, batch_duration: f64, method: ForceCalculationMethod, eta: f64, output_directory: String) -> Result<(), Box<dyn std::error::Error>> {
    clean_output(&output_directory)?;
    set_up_output_directory(&output_directory)?;

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
    //print_all_info(&particles);

    Ok(())
}

fn clean_output(output_directory: &String) -> std::io::Result<()> {
    if Path::new(output_directory).exists() {
        fs::remove_dir_all(output_directory)?;
    }
    Ok(())
}

fn set_up_output_directory(output_directory: &String) -> std::io::Result<()> {
    if !Path::new(output_directory).exists() {
        fs::create_dir_all(output_directory)?;
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