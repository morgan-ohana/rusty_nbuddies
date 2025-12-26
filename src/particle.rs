use serde::{Serialize, Deserialize};

#[derive(PartialEq)]
#[derive(Clone)]
#[derive(Serialize, Deserialize)]
pub struct Particle {
    pub mass: f64,
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub acceleration: [f64;3],
    pub jerk: [f64;3],
    pub snap: [f64;3],
}
