
const SIGMOID_COEFFICIENT: f64 = 4.9;

use petgraph::graphmap::DiGraphMap;
use petgraph::{Graph, Direction};
use petgraph::visit::{self, DfsEvent, EdgeRef};
use std::ops::Range;
use std::collections::HashMap;
use genome::ConnectionGene;

pub struct NeuralNetwork {
    graph: Graph<usize, f64>,
    input_ids: Range<usize>,
}

impl NeuralNetwork {

    pub fn new(inputs: usize, connections: &[ConnectionGene]) -> NeuralNetwork {
        let input_ids = 0 .. inputs;

        // First, find all the nodes we need.
        let mut graph = DiGraphMap::new();
        for conn in connections {
            if conn.enabled {
                graph.add_node(conn.input_node);
                graph.add_node(conn.output_node);
            }
        }

        // Add each connection.
        for conn in connections {
            if conn.enabled {
                graph.add_edge(conn.input_node, conn.output_node, conn.weight);
            }
        }

        // Turn it into a full graph. 
        let full_graph = graph.into_graph();

        NeuralNetwork {
            graph: full_graph,
            input_ids,
        }
    }

    pub fn compute(&self, inputs: &[f64]) -> HashMap<usize, f64> {
        // Using all dangling nodes (nodes with no output connections), we ensure
        // that the depth first search will visit every node in the graph.
        let dangling_nodes: Vec<_> = self.graph.externals(Direction::Outgoing).collect();

        let mut values = HashMap::new();
        for (&value, index) in inputs.iter().zip(self.input_ids.clone()) {
            values.insert(index, value);
        }

        // We're going to first reverse the direction of all the edges, so that
        // nodes point toward the inputs they need. By doing a depth-first search
        // of this, we can ensure that we get a chance to compute the value of a
        // node before using it.
        visit::depth_first_search(&self.graph, dangling_nodes.clone(), |event| {
            match event {
                // When we find a node, give it a dummy value.
                DfsEvent::Discover(id, _) => {
                    values.insert(self.graph[id], 0.0);
                },
                DfsEvent::Finish(id, _) => {
                    // Sum up the values of all the inputs to this node.
                    let total = self.graph.edges_directed(id, Direction::Incoming)
                        .map(|edge| {
                            let source = edge.source();
                            let id = self.graph[source];
                            let value = values[&id];
                            let weight = edge.weight();
                            weight * value
                        })
                        .sum();

                    // Apply the activation function.
                    let value = activation_function(total);

                    // Save this value.
                    values.insert(self.graph[id], value);
                },

                // Ignore edges.
                _ => ()
            }
        });
        
        values
    }
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
        -1.0 * (a.powi(3) * x.powi(3)) / 48.0,
        (a.powi(5) * x.powi(5)) / 480.0,
    ];
    
    terms.into_iter().sum()
}


