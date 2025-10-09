use serde::{Serialize, Deserialize};
use std::fs::{File};
use std::io::{BufWriter, BufReader, Write, Read};
use bincode;
use crate::black_hole::BlackHole;

#[derive(Serialize, Deserialize, Clone)]
pub struct SimulationState {
    pub time: f64,
    pub black_holes: Vec<BlackHole>,
    pub step_count: usize,
    //energy: f64,  // For conservation checking
}

#[derive(Serialize, Deserialize, Clone)]
struct TrajectoryPoint {
    time: f64,
    bh1_x: f64, bh1_y: f64, bh1_vx: f64, bh1_vy: f64,
    bh2_x: f64, bh2_y: f64, bh2_vx: f64, bh2_vy: f64,
}

pub fn save_checkpoint(state: &SimulationState, output_directory: &str, batch_size: &usize) -> anyhow::Result<()> {
    if state.step_count % batch_size != 0 {
	panic!("saving at irregular interval, likely unitended behavior!");
    }
    let file = File::create(format!("{}/restart_{:03}.log", output_directory, state.step_count/batch_size))?;
    let mut writer = BufWriter::new(file);
    
    let encoded = bincode::serialize(state)?;
    writer.write_all(&encoded)?;
    
    println!("Checkpoint saved: (step: {}, time: {:.2})", 
             state.step_count, state.time);
    Ok(())
}

pub fn load_checkpoint(output_directory: &str, batch_num: &usize) -> anyhow::Result<SimulationState> {
    let file = File::open(format!("{}/restart_{:03}.log", output_directory, batch_num))?;
    let mut reader = BufReader::new(file);
    
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;
    
    let state: SimulationState = bincode::deserialize(&buffer)?;
    println!("Checkpoint loaded: (step: {}, time: {:.2})", 
             state.step_count, state.time);
    Ok(state)
}
