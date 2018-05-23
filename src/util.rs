#[cfg(test)]
pub mod test_util {

    use rand::{self, SeedableRng, XorShiftRng};

    pub fn new_rng(maybe_seed: Option<u8>) -> XorShiftRng {
        let seed = maybe_seed.unwrap_or_else(|| rand::random());

        println!(
            "
To reproduce this test, replace:
\tlet mut rng = test_util::new_rng(None)
with:
\tlet mut rng = test_util::new_rng(Some({}))
Don't forget to put the None back when you're done!
",
            seed
        );

        XorShiftRng::from_seed([seed.into(), 1, 1, 1])
    }

}
