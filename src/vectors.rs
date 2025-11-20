
pub fn dot_product(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

pub fn cross_product(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    let mut c: [f64; 3] = [0.0, 0.0, 0.0];
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[2];
    c
}

pub fn magnitude(vec: &[f64; 3]) -> f64 {
    f64::sqrt(dot_product(vec, vec))
}

pub fn add(vec1: &[f64; 3], vec2: &[f64; 3]) -> [f64; 3] {
    let mut sum: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        sum[i] = vec1[i] + vec2[i];
    }
    sum
}

pub fn subtract(vec1: &[f64; 3], vec2: &[f64; 3]) -> [f64; 3] {
    let mut sum: [f64; 3] = [0.0; 3];
    for i in 0..3 {
        sum[i] = vec1[i] - vec2[i];
    }
    sum
}

pub fn scalar_multiply(scalar: &f64, vec: &[f64; 3]) -> [f64; 3] {
    let mut output = vec.clone();
    for i in 0..3 {
        output[i] = scalar * vec[i];
    }
    output
}