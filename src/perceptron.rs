
const SIGMOID_COEFFICIENT: f64 = 4.9;


pub struct Connection {
    pub source: usize,
    pub weight: f64
}

pub struct PerceptronManager<'a> {
    parent: &'a NeuralNetwork,
    values: HashMap<usize, f64>
}

pub struct Perceptron {
    inputs: Vec<Connection>,
}

pub struct NeuralNetwork {
    nodes: Vec<Perceptron>,
}

impl Perceptron {
    pub fn from_inputs(inputs: Vec<Connection>) -> Perceptron {
        Perceptron {
            inputs,
        }
    }

    pub fn get_output(&self, parent: &mut PerceptronManager) -> f64 {
        let total = self.inputs.iter()
            .map(|conn| parent.get_output(conn.source))
            .zip(self.inputs.iter().map(|conn| conn.weight))
            .map(|(output, weight)| output * weight)
            .sum();

        Self::activation_function(total)
    }

    // Uses a taylor series expansion of 1/(1 + e ^ (- a * x))
    // TODO: Use a more efficient approximation algorithm, and/or a lookup table. 
    // We don't need exact values.
    //
    // https://www.wolframalpha.com/input/?i=1%2F(1%2Be%5E(-ax))
    pub fn activation_function(x: f64) -> f64 {
        let a = SIGMOID_COEFFICIENT;
        let terms = vec![
            0.5,
            (a * x) / 4.0,
            -1 * (a.powi(3) * x.powi(3)) / 48.0,
            (a.powi(5) * x.powi(5)) / 480.0,
        ];
        
        terms.sum()
    }
}

impl<'a> PerceptronManager<'a> {
    pub fn new(parent: &NeuralNetwork) -> PerceptronManager {
        PerceptronManager {
            parent,
            values: HashMap::new()
        }
    }

    pub fn 
}

