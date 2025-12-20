use serde::{Serialize, Deserialize};
use std::fs::{File};
use std::io::{BufWriter, BufReader, Write, Read};
use std::path::Path;
use bincode;
use crate::black_hole::BlackHole;

#[derive(Serialize, Deserialize, Clone)]
pub struct SimulationState {
    pub time: f64,
    pub data: Vec<Vec<BlackHole>>,
    pub step_count: usize,
    //energy: f64,  // For conservation checking
}

pub fn save_checkpoint(state: &SimulationState, output_directory: &str, batch_size: &usize) -> anyhow::Result<()> {
    if (state.step_count + 1) % batch_size != 0 {
	    panic!("saving at irregular interval, likely unitended behavior!");
    }
    let batch_num = (state.step_count + 1)/batch_size - 1;
    
    let file = File::create(format!("{}/restart_{:03}.log", output_directory, batch_num))?;
    let mut writer = BufWriter::new(file);        

    let encoded = bincode::serialize(state)?;
    writer.write_all(&encoded)?;

    println!("Checkpoint saved: (num: {}, step: {}, time: {:.2})", 
        batch_num, state.step_count, state.time);
    Ok(())
}

pub fn load_checkpoint(output_directory: &str, batch_num: &usize) -> anyhow::Result<SimulationState> {
    let state: SimulationState;
    
    let file = File::open(format!("{}/restart_{:03}.log", output_directory, batch_num))?;
    let mut reader = BufReader::new(file);
    
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;

    state = bincode::deserialize(&buffer)?;
        
    println!("Checkpoint loaded: (num: {}, step: {}, time: {:.2})", 
        batch_num, state.step_count, state.time);
    Ok(state)
}
