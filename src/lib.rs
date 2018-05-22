extern crate itertools;
extern crate rand;

pub mod genome;
pub mod population;
pub mod util;

use genome::Genome;

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
}

pub struct NeuralNetwork {}
