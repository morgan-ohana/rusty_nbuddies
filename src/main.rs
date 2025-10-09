mod black_hole;
mod forces;
mod time_evol;
mod init_conds;
mod plotting;
mod logging;

use crate::black_hole::BlackHole;
use crate::init_conds::*;
use crate::forces::*;
use crate::time_evol::*;
use crate::plotting::*;
use crate::logging::*;

const AU: f64 = 4.848136811e-9; // AU in kpc

const STEPS: usize = 1000000;
const BATCH_SIZE: usize = 10000;
const DELTA_T: f64 = 0.001*YEAR;
const OUTPUT_DIRECTORY: &str = "output";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Timestep is {}", DELTA_T);

    //let mut black_holes: [BlackHole; 2] = binary_init_conds(10.0*AU, 3.0e-9, 1.0, 2.0);
    let mut black_holes: Vec<BlackHole> = binary_circular_init_conds(1000.0*AU, 10000.0, 10000.0);
    let mut previous_black_holes: Vec<BlackHole> = black_holes.clone();

    let mut data: Vec<Vec<BlackHole>> = vec![black_holes.clone(); STEPS];

    print_all_info(&black_holes);

    for n in 0..STEPS {
        
        data[n] = black_holes.clone();

        recalculate_acceleration_due_to_gravity(&mut black_holes);

        update_velocities(&mut black_holes, &previous_black_holes, &DELTA_T);

        update_positions(&mut black_holes, &previous_black_holes, &DELTA_T);
        
        previous_black_holes = black_holes.clone();

	if n % BATCH_SIZE == 0 {
	    let sim = SimulationState{
		time: n as f64 * DELTA_T,
		black_holes: black_holes.clone(),
		step_count: n
	    };
            save_checkpoint(&sim, &OUTPUT_DIRECTORY, &BATCH_SIZE)?;
	}
    }
    print_all_info(&black_holes);
    create_comprehensive_plots("test", &data, &DELTA_T)?;

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
