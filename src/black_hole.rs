#[derive(PartialEq)]
#[derive(Clone)]
pub struct BlackHole {
    pub mass: f64,
    pub position: [f64; 3],
    pub velocity: [f64; 3],
    pub acceleration: [f64;3],
}