
pub trait NeuralNetwork {
    type Options;

    fn build(genome: &Genome, options: &Self::Options) -> Self;

    fn activate(inputs: &Vec<f64>, output_count: usize) -> Vec<f64>;
}

pub struct PerceptronNetwork {
    nodes: Vec<Perceptron>,
    options: PerceptronOptions
}

pub struct PerceptronOptions {
    max_steps: u64,
}

struct Perceptron {
    inputs: Vec<usize>,
    outputs: Vec<usize>
    previous_value: f64,
    id: usize,
}



