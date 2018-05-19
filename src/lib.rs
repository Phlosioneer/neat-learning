extern crate rand;

mod genome;
use genome::Genome;

pub struct Counter {
    innovation: u64,
    node: usize,
}

impl Counter {
    pub fn new() -> Counter {
        Counter {
            innovation: 0,
            node: 0,
        }
    }

    pub fn next_innovation(&mut self) -> u64 {
        let temp = self.innovation;
        self.innovation += 1;
        temp
    }

    pub fn next_node(&mut self) -> usize {
        let temp = self.node;
        self.node += 1;
        temp
    }
}

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
