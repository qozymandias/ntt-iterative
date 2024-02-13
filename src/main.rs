const MODULUS: u64 = 998_244_353; // Prime modulus for the NTT
const PRIMITIVE_ROOT: u64 = 3; // Primitive root of MODULUS

// Compute (base^exp) % modulus efficiently
fn power_mod(base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result = 1;
    let mut base = base % modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        exp >>= 1;
        base = (base * base) % modulus;
    }
    result
}

fn bit_reverse(idx_in: usize, log2n: usize) -> usize {
    let mut x = idx_in;
    let mut n = 0;
    for i in 0..log2n {
        n <<= 1;
        n |= x & 1;
        x >>= 1;
    }
    n
}

fn ntt(a: &mut [u64], n: usize, primitive_root: u64) {
    let mut m = n;
    let mut h = 0;
    while m > 1 {
        m >>= 1;
        h += 1;
    }
    let mut rev = vec![0; n];
    for i in 0..n {
        rev[i] = rev[i >> 1] >> 1 | (if i & 1 == 1 { n >> 1 } else { 0 });
        if i < rev[i] {
            a.swap(i, rev[i]);
        }
    }
    for i in 1..=h {
        let mh = 1 << i;
        let m = mh >> 1;
        let base = power_mod(primitive_root, (MODULUS - 1) / mh as u64, MODULUS);
        let mut w = 1;
        for j in 0..m {
            for k in (0..n).step_by(mh as usize) {
                let u = a[k + j];
                let t = a[k + j + m] * w % MODULUS;
                a[k + j] = (u + t) % MODULUS;
                a[k + j + m] = (u + MODULUS - t) % MODULUS;
            }
            w = w * base % MODULUS;
        }
    }
}

// Inverse Number Theoretic Transform (NTT)
fn intt(a: &mut [u64], n: usize, primitive_root: u64) {
    let n_inv = power_mod(n as u64, MODULUS - 2, MODULUS);
    ntt(a, n, power_mod(primitive_root, MODULUS - 2, MODULUS));
    for ai in a.iter_mut() {
        *ai = (*ai * n_inv) % MODULUS;
    }
}

fn alt_ntt(a: &mut [u64], n: usize, primitive_root: u64) {
    let mut m = n;
    let mut h = 0;
    while m > 1 {
        m >>= 1;
        h += 1;
    }
    let mut rev = vec![0; n];
    for i in 0..n {
        rev[i] = rev[i >> 1] >> 1 | (if i & 1 == 1 { n >> 1 } else { 0 });
        if i < rev[i] {
            a.swap(i, rev[i]);
        }
    }

    let n = a.len();
    assert!(n.is_power_of_two());
    let mut m = n;
    let mut base = power_mod(primitive_root, (MODULUS - 1) / m as u64, MODULUS);
    while m > 2 {
        m >>= 1;
        let mut r = 0;
        while r < n {
            let mut w = 1;
            for s in r..r + m {
                let u = a[s];
                let d = a[s + m];
                a[s] = (u + d) % MODULUS;
                a[s + m] = (w * ((u + MODULUS - d) % MODULUS)) % MODULUS;
                w = w * base % MODULUS;
            }
            r += 2 * m;
        }
        base = base * base % MODULUS;
    }
    if m > 1 {
        // m = 1
        let mut r = 0;
        while r < n {
            let u = a[r];
            let d = a[r + 1];
            a[r] = (u + d) % MODULUS;
            a[r + 1] = (u + MODULUS - d) % MODULUS;
            r += 2;
        }
    }
}

fn main() {
    {
        let mut coefficients = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let og = coefficients.clone();
        let n = coefficients.len();
        ntt(&mut coefficients, n, PRIMITIVE_ROOT);
        println!("NTT: {:?}", coefficients);

        intt(&mut coefficients, n, PRIMITIVE_ROOT);
        println!("Inverse NTT: {:?}", coefficients);
        assert_eq!(coefficients, og);
    }

    {
        let n = 8;
        let mut vec0: Vec<u64> = vec![4, 1, 4, 2, 1, 3, 5, 6];
        let mut vec1: Vec<u64> = vec![6, 1, 8, 0, 3, 3, 9, 8];
        let expected_out: Vec<u64> = vec![123, 120, 106, 92, 139, 144, 140, 124];

        ntt(&mut vec0, n, PRIMITIVE_ROOT);
        ntt(&mut vec1, n, PRIMITIVE_ROOT);

        let mut res = Vec::with_capacity(n);
        for i in 0..n {
            res.push(vec0[i].clone() * vec1[i].clone());
        }

        intt(&mut res, n, PRIMITIVE_ROOT);
        println!("Inverse NTT: {:?}", res);

        assert_eq!(res, expected_out);
    }

    {
        let mut coefficients = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let og = coefficients.clone();
        let n = coefficients.len();
        alt_ntt(&mut coefficients, n, PRIMITIVE_ROOT);
        println!("NTT: {:?}", coefficients);

        intt(&mut coefficients, n, PRIMITIVE_ROOT);
        println!("Inverse NTT: {:?}", coefficients);
        assert_eq!(coefficients, og);
    }

    {
        let n = 8;
        let mut vec0: Vec<u64> = vec![4, 1, 4, 2, 1, 3, 5, 6];
        let mut vec1: Vec<u64> = vec![6, 1, 8, 0, 3, 3, 9, 8];
        let expected_out: Vec<u64> = vec![123, 120, 106, 92, 139, 144, 140, 124];

        alt_ntt(&mut vec0, n, PRIMITIVE_ROOT);
        alt_ntt(&mut vec1, n, PRIMITIVE_ROOT);

        let mut res = Vec::with_capacity(n);
        for i in 0..n {
            res.push(vec0[i].clone() * vec1[i].clone());
        }

        intt(&mut res, n, PRIMITIVE_ROOT);
        println!("Inverse NTT: {:?}", res);

        assert_eq!(res, expected_out);
    }
}
