use crate::particle::{GravitationalSource, Particle};
use crate::vectors::*;
use crate::forces::{GG, KM_IN_KPC};
use rayon::prelude::*;

#[derive(Debug)]
#[derive(Clone)]
pub enum AccuracyCriterion {
    Geometric(f64), // opening angle theta
    Dynamical(f64), // alpha parameter
}

impl AccuracyCriterion {
    pub fn name(&self) -> String {
        match self {
            AccuracyCriterion::Geometric(theta) => String::from(format!("Geomtric Theta = {theta}")),
            AccuracyCriterion::Dynamical(alpha) => String::from(format!("Dynamical Alpha = {alpha}"))
        }
    }
}

pub enum Node {
    Branch {
        geometric_center: [f64; 3],
        size: f64,
        children: [Option<Box<Node>>; 8],
        mass: f64,
        center_of_mass: [f64; 3],
        velocity_cm: [f64; 3],
        acceleration_cm: [f64; 3],
    },
    Leaf {
        geometric_center: [f64; 3],
        size: f64,
        particle: Particle,
    },
}

impl GravitationalSource for Node {
    fn get_mass(&self) -> f64 {
        match self {
            Node::Branch { mass, .. } => *mass,
            Node::Leaf { particle, .. } => particle.mass,
        }
    }
    fn get_position(&self) -> [f64; 3] {
        match self {
            Node::Branch { center_of_mass, .. } => *center_of_mass,
            Node::Leaf { particle, .. } => particle.position,
        }
    }
    fn get_velocity(&self) -> [f64; 3] {
        match self {
            Node::Branch { velocity_cm, .. } => *velocity_cm,
            Node::Leaf { particle, .. } => particle.velocity,
        }
    }
    fn get_acceleration(&self) -> [f64; 3] {
        match self {
            Node::Branch { acceleration_cm, .. } => *acceleration_cm,
            Node::Leaf { particle, .. } => particle.acceleration,
        }
    }
}

impl Node {
    fn new_branch(geometric_center: [f64; 3], size: f64) -> Self {
        Node::Branch {
            geometric_center,
            size,
            mass: 0.0,
            center_of_mass: [0.0; 3],
            velocity_cm: [0.0; 3],
            acceleration_cm: [0.0; 3],
            children: [None, None, None, None, None, None, None, None],
        }
    }

    fn new_leaf(geometric_center: [f64; 3], size: f64, particle: Particle) -> Self {
        Node::Leaf {
            geometric_center,
            size,
            particle,
        }
    }

    pub fn contains(&self, target: &Particle) -> bool {
        match self {
            Node::Branch { geometric_center, size, .. } | Node::Leaf { geometric_center, size, .. } => {
                for i in 0..3 {
                    if (target.position[i] < geometric_center[i] - size / 2.0) || (target.position[i] > geometric_center[i] + size / 2.0) {
                        return false;
                    }
                }
                true
            }
        }
    }

    pub fn is_approximatable(&self, target: &Particle, previous_target_accel: &[f64; 3], criterion: &AccuracyCriterion) -> bool {
        let distance = {
            let node_pos = self.get_position();
            let mut dist = 0.0;
            for i in 0..3 {
                dist += (node_pos[i] - target.position[i]).powi(2);
            }
            dist.sqrt()
        };
        let size = match self {
            Node::Branch { size, .. } | Node::Leaf { size, .. } => *size,
        };
        match criterion {
            AccuracyCriterion::Geometric(theta) => {
                (size / distance) < *theta
            },
            AccuracyCriterion::Dynamical(alpha) => {
                GG * self.get_mass() * size.powi(2) / distance.powi(4) < alpha * magnitude(previous_target_accel) * KM_IN_KPC
            },
        }
    }

    fn grow_tree(&mut self, particles: Vec<Particle>) {
        match self {
            Node::Leaf { .. } => return, // Leaf nodes do not need to build further
            Node::Branch { geometric_center, size, children, ..} => {            
                // Pre-allocate children
                let mut sorted_particles: [Vec<Particle>; 8] = [Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()];
                let half_size = *size / 2.0;
                
                for particle in particles {
                    let mut octant = 0;
                    for i in 0..3 {
                        if particle.position[i] >= geometric_center[i] {
                            octant |= 1 << i;
                        }
                    }
            
                    sorted_particles[octant].push(particle);
                }

                children.par_iter_mut().zip(sorted_particles).enumerate().for_each(|(octant, (child_opt, particles_for_octant))| {
                    // Do nothing if no particles in this octant
                    if particles_for_octant.is_empty() {
                        return;
                    }

                    // Prepare child center
                    let offset = [
                        if (octant & 1) == 0 { -0.25 } else { 0.25 },
                        if (octant & 2) == 0 { -0.25 } else { 0.25 },
                        if (octant & 4) == 0 { -0.25 } else { 0.25 },
                    ];
                    let child_center = [
                        geometric_center[0] + offset[0] * *size,
                        geometric_center[1] + offset[1] * *size,
                        geometric_center[2] + offset[2] * *size,
                    ];

                    if particles_for_octant.len() == 1 {
                        // Create leaf node if only one particle
                        *child_opt = Some(Box::new(Node::new_leaf(child_center, half_size, particles_for_octant.into_iter().next().unwrap())));
                    } else {
                        // Create branch node
                        let mut child_node = Node::new_branch(child_center, half_size);
                        child_node.grow_tree(particles_for_octant);
                        *child_opt = Some(Box::new(child_node));
                    }
                })
            }
        }
    }

    fn compute_mass_distribution(&mut self) {
        match self {
            Node::Leaf { .. } => return, // Leaf nodes already have mass in their particle
            Node::Branch { geometric_center: _, size: _, mass, center_of_mass, velocity_cm, acceleration_cm, children } => {
                
                // Parallelize work for children
                children.par_iter_mut().for_each(|child_opt| {
                    if let Some(child) = child_opt {
                        child.compute_mass_distribution();
                    }
                });
                
                // Aggregate
                let mut total_mass = 0.0;
                let mut weighted_position = [0.0; 3];
                let mut weighted_velocity = [0.0; 3];
                let mut weighted_acceleration = [0.0; 3];
                
                for child_opt in children.iter_mut() {
                    if let Some(child) = child_opt {
                        let child_mass = child.get_mass();
                        total_mass += child_mass;
                        let child_com = child.get_position();
                        let child_vel_cm = child.get_velocity();
                        let child_accel_cm = child.get_acceleration();

                        for i in 0..3 {
                            weighted_position[i] += child_com[i] * child_mass;
                            weighted_velocity[i] += child_vel_cm[i] * child_mass;
                            weighted_acceleration[i] += child_accel_cm[i] * child_mass;
                        }
                    }
                }

                *mass = total_mass;
                if total_mass > 0.0 {
                    for i in 0..3 {
                        center_of_mass[i] = weighted_position[i] / total_mass;
                        velocity_cm[i] = weighted_velocity[i] / total_mass;
                        acceleration_cm[i] = weighted_acceleration[i] / total_mass;
                    }
                }
            }
        }
    }
}

pub fn build_gravitree(particles: Vec<Particle>) -> Node {
    if particles.is_empty() {
        panic!("Cannot build gravitree with no particles");
    }
    
    let mut root_center: [f64; 3] = [0.0, 0.0, 0.0];
    let mut root_size: f64 = 0.0;
    for particle in &particles {
        for i in 0..3 {
            root_center[i] += particle.position[i];
        }
    }
    for i in 0..3 {
        root_center[i] /= particles.len() as f64;
    }
    for particle in &particles {
        for i in 0..3 {
            root_size = root_size.max((particle.position[i] - root_center[i]).abs() * 2.0);
        }
    }
    root_size *= 1.1; // Slightly enlarge to avoid issues with particles exactly on boundary

    let mut root = Node::new_branch(root_center, root_size);
    root.grow_tree(particles);
    root.compute_mass_distribution();
    root
}