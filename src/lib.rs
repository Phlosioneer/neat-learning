extern crate itertools;
extern crate rand;

pub mod genome;
pub mod population;
pub mod util;

use genome::Genome;
use std::collections::HashMap;

pub struct Counter {
    innovation: u64,
    node: usize,

    /// A mapping from an (input, output) connection to its
    /// corresponding innovation id for this generation.
    ///
    /// This is cleared before and after any evolution steps.
    cached_connections: HashMap<(usize, usize), u64>,

    /// A mapping from an (input, output) being split to its
    /// corresponding (innovation, new_node) ids for this
    /// generation.
    ///
    /// This is cleared before and after any evolution steps.
    cached_splits: HashMap<(usize, usize), (u64, usize)>,
}

impl Counter {
    pub fn new() -> Counter {
        Counter {
            innovation: 0,
            node: 0,
            cached_connections: HashMap::new(),
            cached_splits: HashMap::new(),
        }
    }

    #[cfg(test)]
    pub fn from_id(start: usize) -> Counter {
        Counter {
            innovation: start as u64,
            node: start,
            cached_connections: HashMap::new(),
            cached_splits: HashMap::new(),
        }
    }

    // TODO: Test
    fn next_innovation(&mut self) -> u64 {
        let temp = self.innovation;
        self.innovation += 1;
        temp
    }

    // TODO: Test
    fn next_node(&mut self) -> usize {
        let temp = self.node;
        self.node += 1;
        temp
    }

    /// Returns the innovation number for a new connection between input_node and
    /// output_node. This may be a brand new number, or a reused one from an identical
    /// mutation by another Genome.
    // TODO: Test
    pub fn new_connection(&mut self, input_node: usize, output_node: usize) -> u64 {
        // FIXME: Rust NLL will fix this (no lifetime constraint on the None case).
        {
            let maybe_ret = self.cached_connections.get(&(input_node, output_node));
            if let Some(&innovation) = maybe_ret {
                return innovation;
            }
        }

        let innovation = self.next_innovation();
        self.cached_connections
            .insert((input_node, output_node), innovation);
        innovation
    }

    /// This is called when a genome is trying to split a connection into two
    /// connections with a node between them. It returns the innovation number for
    /// that connection, and the new node's id. These may be brand new, or may be
    /// reused from an identical mutation by another Genome.
    // TODO: Test
    pub fn new_split(&mut self, input_node: usize, output_node: usize) -> (u64, usize) {
        // FIXME: Rust NLL will fix this (no lifetime constraint on the None case).
        {
            let maybe_ret = self.cached_splits.get(&(input_node, output_node));
            if let Some(&(innovation, node)) = maybe_ret {
                return (innovation, node);
            }
        }

        let innovation = self.next_innovation();
        let node = self.next_node();
        self.cached_splits
            .insert((input_node, output_node), (innovation, node));
        (innovation, node)
    }

    // TODO: Test
    pub fn clear_innovations(&mut self) {
        self.cached_connections.clear();
        self.cached_splits.clear();
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
