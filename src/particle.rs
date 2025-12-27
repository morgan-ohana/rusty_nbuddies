use serde::{Serialize, Deserialize};

pub trait GravitationalSource {
    fn get_mass(&self) -> f64;
    fn get_position(&self) -> [f64; 3];
    fn get_velocity(&self) -> [f64; 3];
    fn get_acceleration(&self) -> [f64; 3];
}

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

impl GravitationalSource for Particle {
    fn get_mass(&self) -> f64 {
        self.mass
    }

    fn get_position(&self) -> [f64; 3] {
        self.position
    }

    fn get_velocity(&self) -> [f64; 3] {
        self.velocity
    }

    fn get_acceleration(&self) -> [f64; 3] {
        self.acceleration
    }
}