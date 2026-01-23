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

#[derive(Debug)]
#[derive(Clone)]
pub enum Node {
    Branch {
        bounds: Option<[[f64; 2]; 3]>,
        children: [Box<Node>; 2],
        mass: f64,
        center_of_mass: [f64; 3],
        velocity_cm: [f64; 3],
        acceleration_cm: [f64; 3],
    },
    Leaf {
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

    fn new_branch(children: [Box<Node>; 2]) -> Self {
        Node::Branch {
            bounds: None,
            mass: 0.0,
            center_of_mass: [0.0; 3],
            velocity_cm: [0.0; 3],
            acceleration_cm: [0.0; 3],
            children,
        }
    }

    fn new_leaf(particle: Particle) -> Self {
        Node::Leaf {
            particle
        }
    }

    fn get_bounds(&self) -> [[f64; 2]; 3] {
        match self {
            Node::Leaf { particle } => [
                [particle.position[0], particle.position[0]],
                [particle.position[1], particle.position[1]],
                [particle.position[2], particle.position[2]]
            ],
            Node::Branch { bounds, .. } => bounds.expect("Node::get_bounds should not be called before bounds are initialized!")
        }
    }

    pub fn contains(&self, target: &Particle) -> bool {
        match self {
            Node::Leaf { particle } => target == particle,
            Node::Branch { bounds, ..} => {
                match bounds {
                    None => panic!("Contains check being run on node without initialized bounds"),
                    Some(bounds) => {
                        for i in 0..3 {
                            if target.position[i] < bounds[i][0] || target.position[i] > bounds[i][1] {
                                return false
                            }
                        }
                        return true
                    }
                }
            }
        }
    }

    pub fn is_approximatable(&self, target: &Particle, previous_target_accel: &[f64; 3], criterion: &AccuracyCriterion) -> bool {
        match self {
            Node::Leaf {..} => true,
            Node::Branch {bounds, .. } => {
                let distance = {
                    let node_pos = self.get_position();
                    let mut dist = 0.0;
                    for i in 0..3 {
                        dist += (node_pos[i] - target.position[i]).powi(2);
                    }
                    dist.sqrt()
                };

                //worst case
                let bounds = bounds.expect("Tree is being used for force calculation before it is fully initialized!");
                let size = (bounds[0][1] - bounds[0][0]).min(bounds[1][1] - bounds[1][0]).min(bounds[2][1] - bounds[2][0]);

                match criterion {
                    AccuracyCriterion::Geometric(theta) => {
                        (size / distance) < *theta
                    },
                    AccuracyCriterion::Dynamical(alpha) => {
                        GG * self.get_mass() * size.powi(2) / distance.powi(4) < alpha * magnitude(previous_target_accel) * KM_IN_KPC
                    },
                }
            }
        }
        
        
    }

    fn compute_physical_parameters(&mut self) {
        match self {
            Node::Leaf { .. } => return, // Leaf nodes already have mass in their particle
            Node::Branch { bounds, mass, center_of_mass, velocity_cm, acceleration_cm, children, .. } => {
                
                // Parallelize work for children
                children.par_iter_mut().for_each(|child| {
                    child.compute_physical_parameters();
                });
                
                // Aggregate
                let mut total_mass = 0.0;
                let mut weighted_position = [0.0; 3];
                let mut weighted_velocity = [0.0; 3];
                let mut weighted_acceleration = [0.0; 3];
                
                for child in children.iter_mut() {
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

                // set parameters
                *mass = total_mass;
                *bounds = Some([
                    [children[0].get_bounds()[0][0].min(children[1].get_bounds()[0][0]), children[0].get_bounds()[0][1].max(children[1].get_bounds()[0][1])],
                    [children[0].get_bounds()[1][0].min(children[1].get_bounds()[1][0]), children[0].get_bounds()[1][1].max(children[1].get_bounds()[1][1])],
                    [children[0].get_bounds()[2][0].min(children[1].get_bounds()[2][0]), children[0].get_bounds()[2][1].max(children[1].get_bounds()[2][1])]
                ]);
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

//https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
pub fn build_gravitree(particles: Vec<Particle>) -> Node {
    let (mut leaves, codes) = generate_leaves(particles);

    let branch_structure: Vec<(usize, usize, usize)> = (0..leaves.len() - 1)
        .into_par_iter().map(|idx| {
            let (first, last) = determine_range(&codes, idx);
            let split = find_split(&codes, (first, last));

            (first, last, split)
        }).collect();
    
    let mut root = assemble_tree(leaves.len(), &mut leaves, &branch_structure);

    root.compute_physical_parameters();

    root
}

fn assemble_tree(
    idx: usize,
    leaves: &mut Vec<Option<Node>>,
    branch_structure: &Vec<(usize, usize, usize)>
) -> Node {
    if idx < leaves.len() {
        //this is a leaf, just return itself
        return leaves[idx].take().unwrap()
    }

    // This is an internal Node, it will need it's children before it returns itself
    let node_idx = idx - leaves.len();
    let (first, last, split) = branch_structure[node_idx];

    let left_child_idx = match split == first {
        true => split,
        false => split + leaves.len()
    };

    let right_child_idx = match split + 1 == last {
        true => split + 1,
        false => split + 1 + leaves.len()
    };

    let left_child = assemble_tree(left_child_idx, leaves, branch_structure);
    let right_child = assemble_tree(right_child_idx, leaves, branch_structure);

    Node::new_branch([Box::new(left_child), Box::new(right_child)])
}

//https://developer.nvidia.com/blog/parallelforall/wp-content/uploads/2012/11/karras2012hpg_paper.pdf
fn determine_range(codes: &Vec<u32>, idx: usize) -> (usize, usize) {
    if idx == 0 {
        return (0, codes.len() - 1)
    }
    
    // determine direction or range
    let d = ((codes[idx] ^ codes[idx + 1]).leading_zeros() as isize - (codes[idx] ^ codes[idx - 1]).leading_zeros() as isize).signum();

    let min_prefix = (codes[idx] ^ codes[(idx as isize - d) as usize]).leading_zeros();

    let mut high = match d {
        // Note high should be one larger than the largest possible index such that low can reach the largest possible index if needed
        1 => codes.len() - idx,
        -1 => idx + 1,
        _ => panic!("Direction indicator not 1 or -1. Something wrong with milton codes")
    };
    let mut low = 0;

    while high - low > 1 {
        let mid = (high + low) / 2;

        let mid_prefix = (codes[idx] ^ codes[(idx as isize + (mid as isize)*d) as usize]).leading_zeros();

        if mid_prefix > min_prefix {
            low = mid
        } else {
            high = mid
        }
    }

    match d {
        1 => (idx, idx + low * (d as usize)),
        -1 => ((idx as isize + (low as isize) * d) as usize, idx),
        _ => panic!("Direction indicator not 1 or -1. Something wrong with morton codes")
    }
}

fn find_split(codes: &Vec<u32>, idx_range: (usize, usize)) -> usize {
    let last = idx_range.1;
    let first = idx_range.0;


    //for identical morton codes split in middle
    if codes[first] == codes[last] {
        return (first + last) / 2;
    }

    let common_prefix = (codes[first] ^ codes[last]).leading_zeros();

    let mut high = last;
    let mut low = first;
    while high - low > 1 {
        let mid = (high + low) / 2;

        let mid_prefix = (codes[first] ^ codes[mid]).leading_zeros();

        if mid_prefix > common_prefix {
            low = mid
        } else {
            high = mid
        }
    }

    low
}

fn generate_leaves(particles: Vec<Particle>) -> (Vec<Option<Node>>, Vec<u32>) {
    let bounds = get_bounds(&particles);

    // Generate morton codes
    let mut particle_info: Vec<(Particle, u32)> = particles.into_par_iter().map(|particle| {
        let normalized_position = normalize_coordinates(particle.position, bounds);
        let morton_code = morton_from_normalized(normalized_position);

        (particle, morton_code)
    }).collect();

    // Sort by morton codes
    particle_info.par_sort_by_key(|&(_, code)| code);

    // create leaves
    let (leaves, codes): (Vec<Option<Node>>, Vec<u32>) = particle_info.into_par_iter().map(|(particle, morton_code)| {
        (Some(Node::new_leaf(particle)), morton_code)
    }).collect();

    (leaves, codes)
}

fn get_bounds(particles: &Vec<Particle>) -> [[f64; 2]; 3] {
    let bounds = particles.par_chunks(1024).map(|chunk| {
        let first = &chunk[0];
        let mut local_bounds = [
            [first.position[0], first.position[0]],
            [first.position[1], first.position[1]],
            [first.position[2], first.position[2]]
        ];

        for particle in chunk.iter().skip(1) {
            local_bounds = [
                [local_bounds[0][0].min(particle.position[0]), local_bounds[0][1].max(particle.position[0])],
                [local_bounds[1][0].min(particle.position[1]), local_bounds[1][1].max(particle.position[1])],
                [local_bounds[2][0].min(particle.position[2]), local_bounds[2][1].max(particle.position[2])]
            ];
        }

        local_bounds
    }).reduce(
        || {
            [
                [f64::INFINITY, f64::NEG_INFINITY],
                [f64::INFINITY, f64::NEG_INFINITY],
                [f64::INFINITY, f64::NEG_INFINITY]
            ]
        },
        |bounds_a, bounds_b| {
            [
                [bounds_a[0][0].min(bounds_b[0][0]), bounds_a[0][1].max(bounds_b[0][1])],
                [bounds_a[1][0].min(bounds_b[1][0]), bounds_a[1][1].max(bounds_b[1][1])],
                [bounds_a[2][0].min(bounds_b[2][0]), bounds_a[2][1].max(bounds_b[2][1])]
            ]
        }
    );

    bounds
}

fn expand_bits(v: u32) -> u32 {
    // Spread bits to every 3rd position (inserting two zeros)
    // For 10 bits: 0-9 -> positions 0,3,6,...,27
    let mut x = v as u64;
    x = (x | (x << 16)) & 0x030000FF;
    x = (x | (x << 8)) & 0x0300F00F;
    x = (x | (x << 4)) & 0x030C30C3;
    x = (x | (x << 2)) & 0x09249249;
    x as u32
}


// Calculates a 30-bit Morton code for the given 3D point located within the unit cube [0,1].
fn morton_from_normalized(coords: [f64; 3]) -> u32 {
    // Clamp and scale to [0, 1023]
    let xx = expand_bits((coords[0] * 1024.0).clamp(0.0, 1023.0) as u32);
    let yy = expand_bits((coords[1] * 1024.0).clamp(0.0, 1023.0) as u32);
    let zz = expand_bits((coords[2] * 1024.0).clamp(0.0, 1023.0) as u32);
    
    // Interleave bits: xx * 4 + yy * 2 + zz
    (xx << 2) + (yy << 1) + zz
}

fn normalize_coordinates(coords: [f64; 3], bounds: [[f64; 2]; 3]) -> [f64;3] {
    let mut normalized = [0.0; 3];
    for i in 0..3 {
        normalized[i] = (coords[i] - bounds[i][0]) / (bounds[i][1] - bounds[i][0])
    }
    normalized
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn text_tree_building() {
        let codes = vec![0b00001, 0b00010, 0b00100, 0b00101, 0b10011, 0b11000, 0b11001, 0b11110];

        let branch_structure: Vec<(usize, usize, usize, usize)> = (0..codes.len() - 1)
        .into_iter().map(|idx| {
            let (first, last) = determine_range(&codes, idx);
            let split = find_split(&codes, (first, last));
            (idx, first, last, split)
        }).collect();
        
        let branch_stucture_from_paper = [
            (0,0,7,3),
            (1,0,1,0),
            (2,2,3,2),
            (3,0,3,1),
            (4,4,7,4),
            (5,5,7,6),
            (6,5,6,5)
        ];

        let mut error = 0;
        for i in 0..8 {
            error += (branch_structure[i].0 as isize - branch_stucture_from_paper[i].0 as isize).abs();
            error += (branch_structure[i].1 as isize - branch_stucture_from_paper[i].1 as isize).abs();
            error += (branch_structure[i].2 as isize - branch_stucture_from_paper[i].2 as isize).abs();
            error += (branch_structure[i].3 as isize - branch_stucture_from_paper[i].3 as isize).abs();
        }

        if error != 0 {
            panic!("Something wrong in tree building, likely range finding or split finding")
        }
    }
}