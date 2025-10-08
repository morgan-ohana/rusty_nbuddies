mod black_hole;
mod forces;
mod time_evol;
mod init_conds;
mod plotting;

use crate::black_hole::BlackHole;
use crate::init_conds::*;
use crate::forces::*;
use crate::time_evol::*;
use crate::plotting::*;

const AU: f64 = 4.848136811e-9; // AU in kpc

const STEPS: usize = 100000000;
const DELTA_T: f64 = 0.001*YEAR;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Timestep is {}", DELTA_T);

    //let mut black_holes: [BlackHole; 2] = binary_init_conds(10.0*AU, 3.0e-9, 1.0, 2.0);
    let mut black_holes: [BlackHole; 2] = binary_circular_init_conds(1000.0*AU, 10000.0, 10000.0);
    let mut previous_black_holes: [BlackHole; 2] = black_holes.clone();

    let mut data: Vec<[BlackHole; 2]> = vec![black_holes.clone(); STEPS];

    for n in 0..STEPS {
        
        data[n] = black_holes.clone();

        recalculate_acceleration_due_to_gravity(&mut black_holes);

        update_velocities(&mut black_holes, &previous_black_holes, &DELTA_T);

        update_positions(&mut black_holes, &previous_black_holes, &DELTA_T);

        if n == 0 {
            println!("x pos {}", black_holes[0].position[0]);
            println!("x vel {}", black_holes[0].velocity[0]);
            println!("x accel {}", black_holes[0].acceleration[0]);
            println!("x pos {}", black_holes[1].position[0]);
            println!("x vel {}", black_holes[1].velocity[0]);
            println!("x accel {}", black_holes[1].acceleration[0]);
        }
        
        previous_black_holes = black_holes.clone();
    }

    println!("x pos {}", black_holes[0].position[0]);
    println!("x vel {}", black_holes[0].velocity[0]);
    println!("x accel {}", black_holes[0].acceleration[0]);
    println!("x pos {}", black_holes[1].position[0]);
    println!("x vel {}", black_holes[1].velocity[0]);
    println!("x accel {}", black_holes[1].acceleration[0]);
    create_comprehensive_plots("test", &data, &DELTA_T)?;

    Ok(())
}
