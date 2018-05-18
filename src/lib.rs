extern crate rand;

mod genome;
use genome::Genome;

pub struct Organism {
    genes: Genome,
    brain: NeuralNetwork,
    fitness: f64,
}

pub struct Species {
    members: Vec<Organism>,
}

pub struct Population {
    species: Vec<Species>,
    max_population: usize,
}

pub enum FitnessType {
    Maximize,
    Minimize,
}

pub trait Environment {
    fn fitness_type() -> FitnessType;

    fn calc_fitness(&self, organism: &Organism) -> f64;
}

pub struct NeuralNetwork {}
