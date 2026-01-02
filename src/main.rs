mod particle;
mod driver;
mod gravitree;
mod forces;
mod time_evol;
mod init_conds;
mod plotting;
mod logging;
mod eddington_inverter;
mod vectors;

use std::fs;
use std::path::Path;
use std::env;
use config::Config;

use crate::particle::Particle;
use crate::driver::{run_simulation, GYR};
use crate::gravitree::AccuracyCriterion;
use crate::init_conds::*;
use crate::forces::*;
use crate::plotting::*;
use crate::logging::*;

// Default config values
const DEFAULT_ETA: f64 = 0.01; // timestep accuracy parameter
const DEFAULT_ALPHA: f64 = 1e-3; // dynamical tree accuracy parameter
const DEFAULT_THETA: f64 = 0.3; // geometric tree accuracy parameter 
const DEFAULT_BATCH_DURATION: f64 = 0.1 * GYR; // time per batch
const DEFAULT_MAX_TIME: f64 = 30.0 * GYR;

fn main() {
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
    let init_conds_table_option = settings.get_table("Initial-Conditions");
    
    let init_conds_file_option: Option<String> = match init_conds_table_option {
        Ok(init_conds_table) => Some(match init_conds_table.get("type")
                                            .expect("No type found for initial conditions, ensure you include type = {valid type string} in your initial conditions table")
                                            .clone().into_string().expect("Initial conditions type must be a string, consult documentation for valid types.").as_str() {
            "File" => {
                init_conds_table.get("location")
                    .expect("Unable to read file location, ensure you include location = {some string} in your initial conditions table")
                    .clone().into_string().expect("Location of initial conditions file must be a string!")
            },
            "Binary" => {
                let min_seperation = init_conds_table.get("r_min")
                        .expect("Unable to read minimum seperation, ensure you include r_min = {some float} in your initial conditions table")
                        .clone().into_float().expect("Minimum seperation must be a float!");
                let m1 = init_conds_table.get("m1")
                        .expect("Unable to read particle 1 mass, ensure you include m1 = {some float} in your initial conditions table")
                        .clone().into_float().expect("Particle 1 mass must be a float!");
                let m2 = init_conds_table.get("m2")
                        .expect("Unable to read particle 2 mass, ensure you include m2 = {some float} in your initial conditions table")
                        .clone().into_float().expect("Particle 2 mass must be a float!");
                let eccentricity = match init_conds_table.get("eccentricity") {
                    Some(eccentricity) => eccentricity.clone().into_float().expect("Eccentricity must be a float!"),
                    None => 0.0
                };
                let initial_seperation = match init_conds_table.get("r_init") {
                    Some(initial_seperation) => initial_seperation.clone().into_float().expect("Initial Seperation must be a float!"),
                    None => min_seperation
                };
                let output_path = init_conds_table.get("output-directory")
                        .expect("Ubable to read output path, ensure you include output_directory = {some string} in your inital conditions table")
                        .clone().into_string().expect("Output path must be a string!");
                let particles = binary_init_conds(m1, m2, min_seperation, eccentricity, initial_seperation);
                let file_name = output_path.clone() + "/Binary.log";
                save_init_conds(file_name.clone(), particles).unwrap();
                file_name
            }
            "Plummer" => {
                let r_s = init_conds_table.get("r_s")
                        .expect("Unable to read scale radius, ensure you include r_s = {some float} in your initial conditions table")
                        .clone().into_float().expect("Scale radius must be a float!");
                let total_mass = init_conds_table.get("total_mass")
                        .expect("Unable to read total mass, ensure you include total_mass = {some float} in your initial conditions table")
                        .clone().into_float().expect("Total mass must be a float!");
                let particle_num  = init_conds_table.get("particle_num")
                        .expect("Unable to read particle number, ensure you include particle_num = {some positive int} in your initial conditions table")
                        .clone().into_uint().expect("Particle number must be a positive integer!");
                let output_path = init_conds_table.get("output-directory")
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
                let output_path = init_conds_table.get("output-directory")
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
                let output_path = init_conds_table.get("output-directory")
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
        }),
        _ => None
    };

    // Run the simulation
    if let Ok(simulation_table) = settings.get_table( "Simulation") {
        let init_conds_file = init_conds_file_option.expect("If you wish to run a simulation you must specify initial conditons. Please consult the documentation for how to correclty format the Initial-Condtions table");

        println!("Loading initial conditions from {init_conds_file}");

        let init_conds_data: Vec<Particle> = load_file(init_conds_file).expect("Failed to load initial conditions file, are you sure the path is correct?").data;
        
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
            Some(duration) => duration.clone().into_float().expect("batch duration must be a float!") * GYR,
            None => DEFAULT_BATCH_DURATION
        };

        let max_time = match simulation_table.get("max-time") {
            Some(max_time) => max_time.clone().into_float().expect("Max time must be a float!") * GYR,
            None => DEFAULT_MAX_TIME
        };

        let output_directory = simulation_table.get("output-directory")
            .expect("Ubable to read output path, ensure you include output_directory = {some string} in your simulation table")
            .clone().into_string().expect("Output path must be a string!");

        let mut startup_message = String::from("Begining simulation:");
        startup_message += &(String::from(format!("\n\t output_directory {output_directory}")));
        startup_message += &(String::from(format!("\n\t Calculation Method: ") + &method.name()));
        startup_message += &(String::from(format!("\n\t Timestep Accuracy Parameter {eta}")));
        startup_message += &(String::from(format!("\n\t Batch Duration {}", batch_duration/GYR)));
        startup_message += &(String::from(format!("\n\t Max Time {}", max_time/GYR)));

        println!("{startup_message}");

        clean_output(&output_directory).unwrap();
        set_up_output_directory(&output_directory).unwrap();

        run_simulation(init_conds_data, max_time, batch_duration, method, eta, output_directory);
    }

    // Diagnostics
    if let Ok(plotting_table) = settings.get_table("Diagnostics") {
        let data_directoy = plotting_table.get("data-directory")
            .expect("Must specify data directory for diagnostics to run").clone().into_string().expect("data_directory must be a string");

        let output_directoy = plotting_table.get("output-directory")
            .expect("Must specify output directory for diagnostics to save plots and such too").clone().into_string().expect("output_directory must be a string");

        if let Some(trajectories) = plotting_table.get("trajectories") {
            if trajectories.clone().into_bool().expect("Plotting flag 'trajectories' must be boolean") {plot_trajectories(&data_directoy.clone(), &(output_directoy.clone() + "/trajectories.png")).unwrap()}
        }

        if let Some(energy) = plotting_table.get("energy") {
            if energy.clone().into_bool().expect("Plotting flag 'energy' must be boolean") {plot_energy(&data_directoy.clone(), &(output_directoy.clone() + "/energy.png")).unwrap()}
        }

        if let Some(momentum) = plotting_table.get("momentum") {
            if momentum.clone().into_bool().expect("Plotting flag 'momentum' must be boolean") {plot_momentum(&data_directoy.clone(), &(output_directoy.clone() + "/momentum.png")).unwrap()}
        }

        if let Some(angular_momentum) = plotting_table.get("angular_momentum") {
            if angular_momentum.clone().into_bool().expect("Plotting flag 'angular_momentum' must be boolean") {plot_angular_momentum(&data_directoy.clone(), &(output_directoy.clone() + "/angular_momentum.png")).unwrap()}
        }
    }

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