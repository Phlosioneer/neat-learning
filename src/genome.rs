use itertools::{EitherOrBoth, Itertools};
use rand::Rng;
use std::cmp::Ordering;
use std::slice;

use population::Counter;

type GenomeIter<'a> = slice::Iter<'a, ConnectionGene>;

const CHANCE_TO_INHERIT_DISABLED: f64 = 0.75;
const CHANCE_TO_MUTATE_WEIGHTS: f64 = 0.80;
const CHANCE_TO_ADD_CONNECTION: f64 = 0.50;
const CHANCE_TO_NUDGE_WEIGHT: f64 = 0.50;

#[derive(Clone, Debug, PartialEq)]
pub struct ConnectionGene {
    pub input_node: usize,
    pub output_node: usize,
    pub weight: f64,
    pub enabled: bool,
    pub innovation: u64,
}

impl ConnectionGene {
    pub fn mutate<R: Rng>(&self, mut rng: &mut R, mut counter: &mut Counter, parent: &Genome) -> Vec<ConnectionGene> {
        if !self.enabled {
            return Vec::new()
        }

        if rng.gen_range(0.0, 1.0) < CHANCE_TO_MUTATE_WEIGHTS {
            vec![self.mutate_change_weight(&mut rng)]
        } else {
            if rng.gen_range(0.0, 1.0) < CHANCE_TO_ADD_CONNECTION {
                vec![self.clone(), Self::mutate_add_connection(&mut rng, &mut counter, &parent)]
            } else {
                self.mutate_add_node(&mut counter)
            }
        }
    }

    pub fn mutate_change_weight<R: Rng>(&self, mut rng: &mut R) -> ConnectionGene {
        let mut ret = self.clone();

        if rng.gen_range(0.0, 1.0) < CHANCE_TO_NUDGE_WEIGHT {
            ret.weight += Self::new_weight(&mut rng) / 2.0;
        } else {
            ret.weight = Self::new_weight(&mut rng);
        }

        ret
    }
    
    fn mutate_add_connection<R: Rng>(mut rng: &mut R, counter: &mut Counter, parent: &Genome) -> ConnectionGene {
        let mut nodes = parent.get_nodes();

        let input_node_index = rng.gen_range(0, nodes.len());
        let input_node = nodes.remove(input_node_index);

        let output_node = *rng.choose(&nodes).unwrap();

        let innovation = counter.new_connection(input_node, output_node);

        ConnectionGene {
            input_node,
            output_node,
            weight: Self::new_weight(&mut rng),
            enabled: true,
            innovation,
        }
    }
    
    fn mutate_add_node(&self, counter: &mut Counter) -> Vec<ConnectionGene> {
        let (innovation, new_node) =
            counter.new_split(self.input_node, self.output_node);

        let new_gene = ConnectionGene {
            input_node: new_node,
            output_node: self.output_node,
            weight: 1.0,
            enabled: true,
            innovation,
        };

        let old_gene = ConnectionGene {
            input_node: self.input_node,
            output_node: new_node,
            weight: self.weight,
            enabled: true,
            innovation: self.innovation
        };

        vec![old_gene, new_gene]
    }

    pub fn new_weight<R: Rng>(rng: &mut R) -> f64 {
        rng.gen_range(-1.0, 1.0)
    }

    pub fn crossover<R: Rng>(&self, other: &ConnectionGene, rng: &mut R) -> ConnectionGene {
        let mut ret;
        if rng.gen::<bool>() {
            ret = self.clone();
        } else {
            ret = other.clone();
        }

        if !self.enabled || !other.enabled {
            // 75% chance of being disabled if either parent is disabled.
            // (that includes if they're both disabled! TODO: Remove that;
            // if both parents agree, the child should keep that. Random
            // disabled/enabled flips should be a separate mutation.)
            ret.enabled = !(rng.gen_range(0.0, 1.0) < CHANCE_TO_INHERIT_DISABLED);
        } else {
            ret.enabled = true;
        }

        ret
    }
}

#[derive(Clone, Debug)]
pub struct Genome {
    genes: Vec<ConnectionGene>,
}

impl Genome {
    pub fn new() -> Genome {
        Genome { genes: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Genome {
        Genome {
            genes: Vec::with_capacity(capacity),
        }
    }

    pub fn from_genes(genes: Vec<ConnectionGene>) -> Genome {
        Genome { genes }
    }

    pub fn len(&self) -> usize {
        self.genes.len()
    }

    pub fn add(&mut self, gene: ConnectionGene) {
        self.genes.push(gene);
        self.genes.sort_by_key(|gene| gene.innovation);
    }

    pub fn get_genes(&self) -> &Vec<ConnectionGene> {
        &self.genes
    }

    pub fn iter(&self) -> GenomeIter {
        self.genes.iter()
    }

    pub fn get_nodes(&self) -> Vec<usize> {
        let mut ret = Vec::new();
        for gene in self.genes.iter() {
            if !ret.contains(&gene.input_node) {
                ret.push(gene.input_node);
            }

            if !ret.contains(&gene.output_node) {
                ret.push(gene.output_node);
            }
        }

        ret
    }

    // TODO: Test this.
    pub fn mutate<R: Rng>(&self, mut rng: &mut R, mut counter: &mut Counter) -> Genome {
        let mut copied_genes = self.genes.clone();
        
        let gene_index = rng.gen_range(0, self.genes.len());
        let old_gene: ConnectionGene = copied_genes.remove(gene_index);
        let mutated_genes: Vec<ConnectionGene> = old_gene.mutate(&mut rng, &mut counter, &self);
        
        copied_genes.extend(mutated_genes);
        Genome::from_genes(copied_genes)
    }

    // Do a crossover mutation of the two genomes, where:
    //
    //      fitness(self) <ordering> fitness(other)
    //
    pub fn crossover<R: Rng>(&self, other: &Genome, order: Ordering, mut rng: &mut R) -> Genome {
        let mut ret = Genome::new();

        for gene_match in CrossIter::new(&self, &other) {
            match gene_match {
                GeneMatch::Pair(self_gene, other_gene) => {
                    let gene = self_gene.crossover(other_gene, &mut rng);
                    ret.add(gene);
                }
                GeneMatch::DisjointFirst(gene) => match order {
                    Ordering::Greater | Ordering::Equal => ret.add(gene.clone()),

                    Ordering::Less => (),
                },
                GeneMatch::DisjointSecond(gene) => match order {
                    Ordering::Less | Ordering::Equal => ret.add(gene.clone()),

                    Ordering::Greater => (),
                },
                GeneMatch::ExcessFirst(gene) => match order {
                    Ordering::Greater | Ordering::Equal => ret.add(gene.clone()),

                    Ordering::Less => (),
                },
                GeneMatch::ExcessSecond(gene) => match order {
                    Ordering::Less | Ordering::Equal => ret.add(gene.clone()),

                    Ordering::Greater => (),
                },
            }
        }

        ret
    }

    pub fn new_random<R: Rng>(
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

// An iterator over two genomes that matches up corresponding genes.
#[derive(Debug)]
struct CrossIter<'a> {
    first: GenomeIter<'a>,
    second: GenomeIter<'a>,

    first_peek: Option<&'a ConnectionGene>,
    second_peek: Option<&'a ConnectionGene>,
}

#[derive(Debug)]
enum GeneMatch<'a> {
    Pair(&'a ConnectionGene, &'a ConnectionGene),
    DisjointFirst(&'a ConnectionGene),
    DisjointSecond(&'a ConnectionGene),
    ExcessFirst(&'a ConnectionGene),
    ExcessSecond(&'a ConnectionGene),
}

impl<'a> CrossIter<'a> {
    pub fn new(first: &'a Genome, second: &'a Genome) -> CrossIter<'a> {
        CrossIter {
            first: first.iter(),
            second: second.iter(),

            first_peek: None,
            second_peek: None,
        }
    }
}

impl<'a> Iterator for CrossIter<'a> {
    type Item = GeneMatch<'a>;

    fn next(&mut self) -> Option<GeneMatch<'a>> {
        if self.first_peek.is_none() {
            self.first_peek = self.first.next();
        }

        if self.second_peek.is_none() {
            self.second_peek = self.second.next();
        }

        let first_peek = self.first_peek;
        let second_peek = self.second_peek;
        match (first_peek, second_peek) {
            // No more genes to iterate.
            (None, None) => None,

            // An excess gene.
            (Some(_), None) => Some(GeneMatch::ExcessFirst(self.first_peek.take().unwrap())),
            (None, Some(_)) => Some(GeneMatch::ExcessSecond(self.second_peek.take().unwrap())),

            (Some(better), Some(worse)) => match better.innovation.cmp(&worse.innovation) {
                // If they're equal, they're a pair.
                Ordering::Equal => Some(GeneMatch::Pair(
                    self.first_peek.take().unwrap(),
                    self.second_peek.take().unwrap(),
                )),

                // A disjoint gene. Return the one with the lower innovation number.
                Ordering::Less => Some(GeneMatch::DisjointFirst(self.first_peek.take().unwrap())),
                Ordering::Greater => {
                    Some(GeneMatch::DisjointSecond(self.second_peek.take().unwrap()))
                }
            },
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::XorShiftRng;
    use std::cmp;
    use util::test_util;

    #[test]
    fn test_from_genes() {
        let (test1, test2, test3) = paper_crossover_example();

        do_test_from_genes(test1);
        do_test_from_genes(test2);
        do_test_from_genes(test3);
    }

    fn do_test_from_genes(test_genes: Vec<ConnectionGene>) {
        let genome = Genome::from_genes(test_genes.clone());
        assert_eq!(genome.get_genes(), &test_genes);
    }

    #[test]
    fn test_add_gene() {
        let (test1, test2, test3) = paper_crossover_example();

        do_test_add_gene(test1);
        do_test_add_gene(test2);
        do_test_add_gene(test3);
    }

    fn do_test_add_gene(test_genes: Vec<ConnectionGene>) {
        let mut genome = Genome::new();

        for gene in test_genes.iter() {
            genome.add(gene.clone());
        }

        assert_eq!(genome.get_genes(), &test_genes);
    }

    #[test]
    fn test_genome_len() {
        let (test1, test2, test3) = paper_crossover_example();

        do_test_genome_len(test1);
        do_test_genome_len(test2);
        do_test_genome_len(test3);
    }

    fn do_test_genome_len(test_genes: Vec<ConnectionGene>) {
        let genome = Genome::from_genes(test_genes.clone());
        assert_eq!(genome.len(), test_genes.len());
    }

    #[test]
    fn test_get_nodes() {
        let (genes, _, _) = paper_crossover_example();

        let genome = Genome::from_genes(genes);
        let nodes = genome.get_nodes();

        for i in 1..6 {
            assert!(nodes.contains(&i), "Node id not found: {}", i);
        }

        assert!(nodes.contains(&8), "Node id not found: 8");

        assert_eq!(nodes.len(), 6);
    }

    #[test]
    fn test_add_node() {
        let (test1, test2, test3) = paper_crossover_example();

        do_test_add_node(test1);
        do_test_add_node(test2);
        do_test_add_node(test3);
    }

    fn do_test_add_node(test_genes: Vec<ConnectionGene>) {
        let parent = Genome::from_genes(test_genes.clone());

        for gene in test_genes {
            println!("Testing gene: {:#?}", gene);
            let mut counter = Counter::from_id(100);
            let output_genes = gene.mutate_add_node(&mut counter);
            println!("Output: {:#?}\n", output_genes);

            assert_eq!(output_genes.len(), 2, "Wrong number of output genes.");

            let out1 = output_genes[0].clone();
            let out2 = output_genes[1].clone();

            let start;
            let end;
            if out1.output_node == out2.input_node {
                start = out1;
                end = out2;
            } else if out2.output_node == out1.input_node {
                start = out2;
                end = out1;
            } else {
                panic!("The genes don't have any inputs/outputs in common");
            }

            assert_eq!(start.input_node, gene.input_node, "Input node doesn't match.");
            assert_eq!(end.output_node, gene.output_node, "Output node doesn't match.");
            assert_eq!(start.enabled && end.enabled, true, "A node is disabled.");
            assert_eq!(start.output_node, 100, "The new node must have a new id.");
            if start.weight == 1.0 {
                assert_eq!(end.weight, gene.weight, "Original weight was not preserved.");
                assert_eq!(end.innovation, gene.innovation, "Original innovation was not preserved or is backwards.");
                assert_eq!(start.innovation, 100, "One of the nodes must be a new innovation number.");
            } else if end.weight == 1.0 {
                assert_eq!(start.weight, gene.weight, "Original weight was not preserved.");
                assert_eq!(start.innovation, gene.innovation, "Original innovation was not preserved or is backwards.");
                assert_eq!(end.innovation, 100, "One of the nodes must be a new innovation number.");
            } else {
                panic!("At least one connection must have a weight of 1.0");
            }
        }
    }

    #[test]
    fn test_add_connection() {
        let (test1, test2, test3) = paper_crossover_example();

        let mut rng = test_util::new_rng(Some(212));

        do_test_add_connection(test1, &mut rng);
        do_test_add_connection(test2, &mut rng);
        do_test_add_connection(test3, &mut rng);
    }

    fn do_test_add_connection(
        test_genome: Vec<ConnectionGene>,
        mut rng: &mut XorShiftRng,
    ) {
        let parent = Genome::from_genes(test_genome);
        let nodes = parent.get_nodes();

        for _ in 0 .. 100 {
            let mut counter = Counter::from_id(100);
            let output_gene = ConnectionGene::mutate_add_connection(&mut rng, &mut counter, &parent);
            println!("output gene: {:#?}\n", output_gene);
            
            assert!(nodes.contains(&output_gene.input_node), "Input node ID doesn't exist. Node ids: {:?}", nodes);
            assert!(nodes.contains(&output_gene.output_node), "Output node ID doesn't exist. Node ids: {:?}", nodes);
            
            assert_eq!(output_gene.enabled, true, "Output gene is disabled.");
            assert_eq!(output_gene.innovation, 100, "Output gene has incorrect innovation number (expected 100)");
        }
    }

    #[test]
    fn test_crossover() {
        let (in1, in2, out) = paper_crossover_example();

        do_test_crossover(in1.clone(), in2.clone(), out.clone(), Ordering::Equal);
        do_test_crossover(in2.clone(), in1.clone(), out.clone(), Ordering::Equal);

        do_test_crossover(in1.clone(), in2.clone(), in1.clone(), Ordering::Greater);
        do_test_crossover(in2.clone(), in1.clone(), in2.clone(), Ordering::Greater);

        do_test_crossover(in1.clone(), in2.clone(), in2.clone(), Ordering::Less);
        do_test_crossover(in2.clone(), in1.clone(), in1.clone(), Ordering::Less);
    }

    fn do_test_crossover(
        in1: Vec<ConnectionGene>,
        in2: Vec<ConnectionGene>,
        out: Vec<ConnectionGene>,
        order: Ordering,
    ) {
        println!(
            "input: {:#?}\n{:#?}\noutput: {:#?}\norder: {:?}",
            &in1, &in2, &out, order
        );

        let mut rng = test_util::new_rng(None);

        let genome1 = Genome::from_genes(in1.clone());
        let genome2 = Genome::from_genes(in2.clone());
        let output_genome = genome1.crossover(&genome2, order, &mut rng);

        assert_eq!(output_genome.get_genes(), &out);
    }

    #[test]
    fn test_random_choice_of_enabled() {
        let in1 = vec![ConnectionGene {
            input_node: 1,
            output_node: 2,
            weight: 0.0,
            enabled: false,
            innovation: 1,
        }];

        let in2 = vec![ConnectionGene {
            input_node: 1,
            output_node: 2,
            weight: 0.0,
            enabled: true,
            innovation: 1,
        }];

        let mut rng = test_util::new_rng(None);

        do_test_random_choice_of_enabled(&in1, &in2, &mut rng);
        do_test_random_choice_of_enabled(&in2, &in1, &mut rng);

        // In the original paper, it could be enabled even if both are disabled.
        do_test_random_choice_of_enabled(&in1, &in1, &mut rng);
    }

    fn do_test_random_choice_of_enabled(
        in1: &Vec<ConnectionGene>,
        in2: &Vec<ConnectionGene>,
        mut rng: &mut XorShiftRng,
    ) {
        let mut false_output_found = false;
        let mut true_output_found = false;
        for _ in 0..100 {
            let genome1 = Genome::from_genes(in1.clone());
            let genome2 = Genome::from_genes(in2.clone());
            let output = genome1.crossover(&genome2, Ordering::Greater, &mut rng);
            if output.get_genes()[0].enabled {
                true_output_found = true;
            } else {
                false_output_found = true;
            }

            if true_output_found && false_output_found {
                break;
            }
        }

        assert!(
            true_output_found && false_output_found,
            "true found: {}\nfalse_found: {}",
            true_output_found,
            false_output_found
        );
    }

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
                let genome = Genome::new_random(&mut rng, &mut counter, inputs, outputs);

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
                let genome = Genome::new_random(&mut rng, &mut counter, inputs, outputs);

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

    // This is the same example the paper uses on p109, but gene #5 is not disabled
    // in either parent. This method isn't testing the randomness.
    fn paper_crossover_example() -> (
        Vec<ConnectionGene>,
        Vec<ConnectionGene>,
        Vec<ConnectionGene>,
    ) {
        let first = vec![
            ConnectionGene {
                input_node: 1,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 1,
            },
            ConnectionGene {
                input_node: 2,
                output_node: 4,
                weight: 0.0,
                enabled: false,
                innovation: 2,
            },
            ConnectionGene {
                input_node: 3,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 3,
            },
            ConnectionGene {
                input_node: 2,
                output_node: 5,
                weight: 0.0,
                enabled: true,
                innovation: 4,
            },
            ConnectionGene {
                input_node: 5,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 5,
            },
            ConnectionGene {
                input_node: 1,
                output_node: 8,
                weight: 0.0,
                enabled: true,
                innovation: 8,
            },
        ];

        let second = vec![
            ConnectionGene {
                input_node: 1,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 1,
            },
            ConnectionGene {
                input_node: 2,
                output_node: 4,
                weight: 0.0,
                enabled: false,
                innovation: 2,
            },
            ConnectionGene {
                input_node: 3,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 3,
            },
            ConnectionGene {
                input_node: 2,
                output_node: 5,
                weight: 0.0,
                enabled: true,
                innovation: 4,
            },
            ConnectionGene {
                input_node: 5,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 5,
            },
            ConnectionGene {
                input_node: 5,
                output_node: 6,
                weight: 0.0,
                enabled: true,
                innovation: 6,
            },
            ConnectionGene {
                input_node: 6,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 7,
            },
            ConnectionGene {
                input_node: 3,
                output_node: 5,
                weight: 0.0,
                enabled: true,
                innovation: 9,
            },
            ConnectionGene {
                input_node: 1,
                output_node: 6,
                weight: 0.0,
                enabled: true,
                innovation: 10,
            },
        ];

        let third = vec![
            ConnectionGene {
                input_node: 1,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 1,
            },
            ConnectionGene {
                input_node: 2,
                output_node: 4,
                weight: 0.0,
                enabled: false,
                innovation: 2,
            },
            ConnectionGene {
                input_node: 3,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 3,
            },
            ConnectionGene {
                input_node: 2,
                output_node: 5,
                weight: 0.0,
                enabled: true,
                innovation: 4,
            },
            ConnectionGene {
                input_node: 5,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 5,
            },
            ConnectionGene {
                input_node: 5,
                output_node: 6,
                weight: 0.0,
                enabled: true,
                innovation: 6,
            },
            ConnectionGene {
                input_node: 6,
                output_node: 4,
                weight: 0.0,
                enabled: true,
                innovation: 7,
            },
            ConnectionGene {
                input_node: 1,
                output_node: 8,
                weight: 0.0,
                enabled: true,
                innovation: 8,
            },
            ConnectionGene {
                input_node: 3,
                output_node: 5,
                weight: 0.0,
                enabled: true,
                innovation: 9,
            },
            ConnectionGene {
                input_node: 1,
                output_node: 6,
                weight: 0.0,
                enabled: true,
                innovation: 10,
            },
        ];

        (first, second, third)
    }

}
