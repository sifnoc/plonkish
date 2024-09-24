use halo2_curves::{bn256::G1Affine, pairing::Engine, serde::SerdeObject};
use itertools::Itertools;
use plonkish_backend::{
    halo2_curves::{bn256::Bn256, group::ff::Field},
    util::arithmetic::{
        batch_projective_to_affine, fixed_base_msm, powers, radix2_fft, root_of_unity_inv,
        window_size, window_table, PrimeField,
    },
};
use rand::rngs::OsRng;
use std::{env, fs::File, io::Write};

// Some of code and logic are referenced from `https://github.com/han0110/halo2-kzg-srs`
fn main() {
    let dst_prefix = env::args()
        .nth(1)
        .expect("Please specify destination file path prefix (will be appended with suffix k)");
    let desired_k = env::args()
        .nth(2)
        .and_then(|s| s.parse::<u32>().ok())
        .expect("Please specify the number of K");

    // Generate destination file
    //
    // The logic is referenced from the `src/pcs/univariate/kzg.rs` file
    let num_vars = desired_k as usize;
    let poly_size = 1 << num_vars;
    let ws = window_size(poly_size);

    let g1 = <Bn256 as Engine>::G1Affine::generator();
    let g2 = <Bn256 as Engine>::G2Affine::generator();

    let window_table_g1 = window_table::<G1Affine>(ws, g1);
    let s = <Bn256 as Engine>::Scalar::random(OsRng);

    let monomial = powers(s).take(poly_size).collect_vec();

    let monomial_g1: Vec<<Bn256 as Engine>::G1Affine> =
        batch_projective_to_affine(&fixed_base_msm(ws, &window_table_g1, &monomial));

    let lagrange_g1: Vec<<Bn256 as Engine>::G1Affine> = {
        let n_inv = <Bn256 as Engine>::Scalar::TWO_INV.pow_vartime([num_vars as u64]);
        let mut lagrange = monomial;
        radix2_fft(&mut lagrange, root_of_unity_inv(num_vars), num_vars);
        lagrange.iter_mut().for_each(|v| *v *= n_inv);
        batch_projective_to_affine(&fixed_base_msm(ws, &window_table_g1, &lagrange))
    };

    let powers_of_s_g2: Vec<<Bn256 as Engine>::G2Affine> = {
        let powers_of_s_g2 = powers(s).take(poly_size).collect_vec();
        let window_table_g2 = window_table(ws, g2);
        batch_projective_to_affine(&fixed_base_msm(ws, &window_table_g2, &powers_of_s_g2))
    };

    // Exports the SRS to a file
    let mut writer = File::create(format!("{}{}", dst_prefix, num_vars)).unwrap();

    writer.write_all(&desired_k.to_le_bytes()).unwrap();

    // Writes them to the file
    for m in monomial_g1.iter() {
        m.write_raw(&mut writer).unwrap();
    }
    for l in lagrange_g1.iter() {
        l.write_raw(&mut writer).unwrap();
    }
    for p in powers_of_s_g2.iter() {
        p.write_raw(&mut writer).unwrap();
    }

    println!("SRS generated successfully");
}
