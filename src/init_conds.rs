use crate::black_hole::BlackHole;
use crate::forces::GG;
use crate::forces::KM_PER_KPC;

pub fn binary_init_conds(seperation: f64, angular_momentum: f64, m1: f64, m2: f64) -> [BlackHole; 2] {
    // 0 = p1 + p2
    // l = p1 r1 + p2 r2 = p1 (r1 - r2) = p1 sep 
    // v1 = l/m1 sep
    // v2 = -l/m2 sep
    let bh_1 = BlackHole {
        mass: m1,
        position: [0.5*seperation, 0.0, 0.0],
        velocity: [0.0, angular_momentum / (m1 * seperation), 0.0],
        acceleration: [0.0, 0.0, 0.0],
    };
    let bh_2 = BlackHole {
        mass: m2,
        position: [-0.5*seperation, 0.0, 0.0],
        velocity: [0.0, - angular_momentum / (m2 * seperation), 0.0],
        acceleration: [0.0, 0.0, 0.0],
    };
    [bh_1, bh_2]
}

pub fn binary_circular_init_conds(seperation: f64, m1: f64, m2: f64) -> [BlackHole; 2] {
    // 0 = m1*x + m2 (x - sep) = - m2 sep + (m1 + m2)x => x = m2 sep/ (m1 + m2)
    let offset = m2 * seperation / (m1+ m2);
    // G (m1m2/mu) / sep^2 = G (m1 + m2) / sep ^2= a = omega^2 sep => omega^2 = G mu / sep^3 => v_i = r_i omega = r_i sqrt(G (m1+ m2) / sep^3); [r_i G (m1+m2) / sep^3]) = kpc (km^2 kpc Msun / (Msun s^2 kpc^3))^1/2 = km/s
    // [omega] = [sqrt(G (m1+ m2) / sep^3)] = [km^2 kpc Msun / Msun s^2 kpc^3]^1/2 = km / s kpc
    let omega = f64::sqrt(GG*(m1+m2)/(seperation*seperation*seperation))*1.3; //1.1 fudge factor to make it a tad elliptical
    println!("ensure timestep less than 2/omega = {}", (2.0*KM_PER_KPC/omega));

    let bh_1 = BlackHole {
        mass: m1,
        position: [offset, 0.0, 0.0],
        velocity: [0.0, offset * omega, 0.0],
        acceleration: [0.0, 0.0, 0.0],
    };
    let bh_2 = BlackHole {
        mass: m2,
        position: [offset - seperation, 0.0, 0.0],
        velocity: [0.0, (offset - seperation) * omega, 0.0],
        acceleration: [0.0, 0.0, 0.0],
    };
    [bh_1, bh_2]
}


