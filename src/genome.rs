use rand::Rng;
use std::cmp::{self, Ordering};
use std::slice;

type GenomeIter<'a> = slice::Iter<'a, ConnectionGene>;

const CHANCE_TO_INHERIT_DISABLED: f64 = 0.75;

#[derive(Clone, Debug, PartialEq)]
pub struct ConnectionGene {
    pub input_node: usize,
    pub output_node: usize,
    pub weight: f64,
    pub enabled: bool,
    pub innovation: u64,
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

    // Do a crossover mutation of the two genomes, where:
    //
    //      fitness(self) <ordering> fitness(other)
    //
    pub fn crossover<R: Rng>(&self, other: &Genome, order: Ordering, rng: &mut R) -> Genome {
        let max_len = cmp::max(self.len(), other.len());

        let mut ret = Genome::new();

        for gene_match in CrossIter::new(&self, &other) {
            println!("Match: {:?} Ordering: {:?}", gene_match, order);
            match gene_match {
                GeneMatch::Pair(self_gene, other_gene) => {
                    let mut gene;
                    if rng.gen::<bool>() {
                        gene = self_gene.clone();
                    } else {
                        gene = other_gene.clone();
                    }

                    if !self_gene.enabled || !other_gene.enabled {
                        // 75% chance of being disabled if either parent is disabled.
                        // (that includes if they're both disabled! TODO: Remove that;
                        // if both parents agree, the child should keep that. Random
                        // disabled/enabled flips should be a separate mutation.)
                        gene.enabled = !(rng.gen_range(0.0, 1.0) < CHANCE_TO_INHERIT_DISABLED);
                    } else {
                        gene.enabled = true;
                    }

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
    use rand::{SeedableRng, XorShiftRng};

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

        let genome1 = Genome::from_genes(in1.clone());
        let genome2 = Genome::from_genes(in2.clone());
        let mut dummy_rand = XorShiftRng::from_seed([1, 1, 1, 1]);
        let output_genome = genome1.crossover(&genome2, order, &mut dummy_rand);

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

        let mut rng = XorShiftRng::from_seed([1, 1, 1, 1]);

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
