use genome::{ConnectionGene, Genome};
use rand::Rng;
use Counter;
use Environment;

pub struct Population<R: Rng> {
    members: Vec<Genome>,
    max_size: usize,
    rng: R,
    counter: Counter,
}

impl<R: Rng> Population<R> {
    pub fn new(max_size: usize, rng: R) -> Population<R> {
        Population {
            members: Vec::new(),
            max_size,
            rng,
            counter: Counter::new(),
        }
    }

    /// Creates the initial population. If there is already a population, it
    /// deletes it and regenerates a new one.
    // TODO: Test that the population is diverse
    pub fn initialize<E: Environment>(&mut self, env: &mut E) {
        self.members.clear();

        let input_count = env.input_count();
        let output_count = env.output_count();

        if input_count == 0 {
            panic!("0 inputs to genome! Check your Environment impl.");
        }

        if output_count == 0 {
            panic!("0 outputs from genome! Check your Environment impl.");
        }

        self.counter.clear_innovations();
        for _ in 0..self.max_size {
            let member =
                Genome::new_random(&mut self.rng, &mut self.counter, input_count, output_count);
            self.members.push(member);
        }
        self.counter.clear_innovations();
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn max_len(&self) -> usize {
        self.max_size
    }

    /*pub fn evolve_once(&mut self, env: &mut Environment) -> Genome {
       
    }*/
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::XorShiftRng;
    use util::test_util;
    use {FitnessType, Organism};

    struct TestEnv;

    impl Environment for TestEnv {
        fn input_count(&self) -> usize {
            10
        }

        fn output_count(&self) -> usize {
            4
        }

        fn fitness_type() -> FitnessType {
            panic!()
        }

        fn calc_fitness(&self, _: &Organism) -> f64 {
            panic!()
        }
    }

    #[test]
    pub fn test_initialize_makes_max_size() {
        let mut rng = test_util::new_rng(None);

        let mut pop = Population::new(150, rng);
        pop.initialize(&mut TestEnv);

        assert_eq!(pop.len(), 150);
    }
}
