extern crate itertools;
extern crate rand;
extern crate petgraph;

pub mod genome;
pub mod population;
pub mod util;
pub mod perceptron;

use genome::{ConnectionGene, Genome};
use itertools::{EitherOrBoth, Itertools};
use population::Counter;
use rand::Rng;
use perceptron::NeuralNetwork;

pub struct Organism {
    genes: Genome,
    brain: NeuralNetwork,
    fitness: f64,
}

pub struct Species {
    members: Vec<Organism>,
}

pub enum FitnessType {
    Maximize,
    Minimize,
}

pub trait Environment {
    fn fitness_type() -> FitnessType;

    fn input_count(&self) -> usize;
    fn output_count(&self) -> usize;

    fn calc_fitness(&self, organism: &Organism) -> f64;
    fn random_genome<R: Rng>(
        rng: &mut R,
        counter: &mut Counter,
        input_count: usize,
        output_count: usize,
    ) -> Genome;
}

pub struct NetworkEnv;

impl Environment for NetworkEnv {
    fn fitness_type() -> FitnessType {
        panic!()
    }

    fn input_count(&self) -> usize {
        10
    }

    fn output_count(&self) -> usize {
        4
    }

    fn calc_fitness(&self, _: &Organism) -> f64 {
        panic!()
    }

    fn random_genome<R: Rng>(
        mut rng: &mut R,
        counter: &mut Counter,
        input_count: usize,
        output_count: usize,
    ) -> Genome {
        // We need to supply at least enough connections to define each input
        // and output node.
        //
        // This code sets up a pair of iterators over the possible input node
        // id's and output node id's in a random order. Then, it zips them together,
        // and uses the zipped pairs as connection values. This will uniquely connect
        // an input and an output.
        //
        // Any excess (because the inputs and outputs are a different size) are
        // chosen randomly from the appropriate input/output ids.
        //
        // This process ensures that each input and output have been used at least
        // once, and that the minimum number of connections are made.
        let mut inputs: Vec<usize> = (0..input_count).collect();
        let mut outputs: Vec<usize> = (0..output_count).map(|x| x + input_count).collect();

        let remainder = if inputs > outputs {
            inputs.clone()
        } else {
            outputs.clone()
        };

        rng.shuffle(&mut inputs);
        rng.shuffle(&mut outputs);

        let mut genes = Vec::new();
        for zip in inputs.into_iter().zip_longest(outputs.into_iter()) {
            let (input, output) = match zip {
                EitherOrBoth::Both(input, output) => (input, output),
                EitherOrBoth::Left(input) => (input, *rng.choose(&remainder).unwrap()),
                EitherOrBoth::Right(output) => (*rng.choose(&remainder).unwrap(), output),
            };

            // If two genomes would make a new connection between the same two
            // nodes during an evolution step, they must share the same innovation
            // number.
            let innovation = counter.new_connection(input, output);

            let gene = ConnectionGene {
                input_node: input,
                output_node: output,
                weight: ConnectionGene::new_weight(&mut rng),
                enabled: true,
                innovation,
            };

            genes.push(gene);
        }

        Genome::from_genes(genes)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::cmp;
    use util::test_util;

    // This ensures the random genomes have the absolute minimum number of connections
    // required to define every input/output node id.
    #[test]
    fn test_new_random_minimum_length() {
        let mut rng = test_util::new_rng(None);

        for inputs in 1..10 {
            for outputs in 1..10 {
                println!("Testing {} inputs and {} outputs.", inputs, outputs);
                let minimum_length = cmp::max(inputs, outputs);

                let mut counter = Counter::new();
                let genome = NetworkEnv::random_genome(&mut rng, &mut counter, inputs, outputs);

                assert_eq!(
                    genome.len(),
                    minimum_length,
                    "Genome is not minimal: {:#?}",
                    genome
                );
            }
        }
    }

    // This ensures the random genomes define every input/output node id.
    #[test]
    fn test_new_random_node_ids() {
        let mut rng = test_util::new_rng(None);

        for inputs in 1..10 {
            for outputs in 1..10 {
                println!("Testing {} inputs and {} outputs.", inputs, outputs);

                let mut counter = Counter::new();
                let genome = NetworkEnv::random_genome(&mut rng, &mut counter, inputs, outputs);

                let nodes = genome.get_nodes();

                for id in 0..(inputs + outputs) {
                    assert!(
                        nodes.contains(&id),
                        "Id {} not found in genome: {:#?}",
                        id,
                        genome
                    );
                }
            }
        }
    }
}
