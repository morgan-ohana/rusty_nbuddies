use serde::{Serialize, Deserialize};
use std::fs::{File};
use std::io::{BufWriter, BufReader, Write, Read};
use bincode;
use crate::particle::Particle;

#[derive(Serialize, Deserialize, Clone)]
pub struct SimulationState {
    pub time: f64,
    pub data: Vec<Particle>,
    pub step_count: usize,
}

pub fn save_checkpoint(state: &SimulationState, output_directory: &str, batch_num: &usize) -> anyhow::Result<()> {
    let file = File::create(format!("{}/snapshot_{:03}.log", output_directory, batch_num))?;
    let mut writer = BufWriter::new(file);        

    let encoded = bincode::serialize(state)?;
    writer.write_all(&encoded)?;

    println!("Checkpoint saved: (num: {}, step: {}, time: {:.2})", 
        batch_num, state.step_count, state.time);
    Ok(())
}

pub fn load_checkpoint(output_directory: &str, batch_num: &usize) -> anyhow::Result<SimulationState> {
    let state: SimulationState;
    
    let file = File::open(format!("{}/snapshot_{:03}.log", output_directory, batch_num))?;
    let mut reader = BufReader::new(file);
    
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;

    state = bincode::deserialize(&buffer)?;
        
    println!("Checkpoint loaded: (num: {}, step: {}, time: {:.2})", 
        batch_num, state.step_count, state.time);
    Ok(state)
}

pub fn load_file(file_name: String) -> anyhow::Result<SimulationState> {
    let state: SimulationState;
    
    let file = File::open(format!("{file_name}"))?;
    let mut reader = BufReader::new(file);
    
    let mut buffer = Vec::new();
    reader.read_to_end(&mut buffer)?;

    state = bincode::deserialize(&buffer)?;
        
    println!("File loaded from {file_name}");
    Ok(state)
}

pub fn save_init_conds(file_name: String, data: Vec<Particle>) -> anyhow::Result<()> {
    let state: SimulationState = SimulationState { time: 0.0, data: data , step_count: 0 };

    let file = File::create(format!("{file_name}"))?;
    let mut writer = BufWriter::new(file);        

    let encoded = bincode::serialize(&state)?;
    writer.write_all(&encoded)?;

    println!("Saved Initial Conditions at {file_name}");
    Ok(())
}
