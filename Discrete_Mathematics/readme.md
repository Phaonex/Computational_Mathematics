# Computational Mathematics: Discrete Mathematics

This repository contains rigorous computational implementations of Discrete Mathematics, focusing on the structures that form the basis of advanced data science, numerical analysis, and machine learning[cite: 1]. 

All logic is implemented exclusively in Julia, which serves as the primary high-performance computing (HPC) language for these mathematical models, positioned to replace Node.js or FastAPI wherever backend bottlenecks occur.

## Top 10 Core Lectures

The following implementations represent the most critical milestones in the curriculum, focusing heavily on combinatorics, algorithm complexity, cryptography, and network topology[cite: 1]:

*   **Lecture 2: Basic Concepts of Combinatorics** (`lecture_2_basic_concepts_of_combinatorics.jl`)
    *   *Focus:* Fundamental subset, sequence, and permutation architectures underpinning probability spaces[cite: 1].
*   **Lecture 6: The Principle of Inclusion-Exclusion** (`lecture_6_the_principle_of_inclusion_exclusion.jl`)
    *   *Focus:* Advanced counting techniques and Stirling numbers for complex set intersections[cite: 1].
*   **Lecture 8: Linear Recurrences and Fibonacci Numbers** (`lecture_8_linear_recurrences_and_fibonacci_numbers.jl`)
    *   *Focus:* Polynomial methods for solving recurrent algorithmic sequences[cite: 1].
*   **Lecture 10: The Structure of Numbers** (`lecture_10_the_structure_of_numbers.jl`)
    *   *Focus:* Binary representations and the fundamental building blocks of computational memory[cite: 1].
*   **Lecture 13: Enormous Exponents and Card Shuffling** (`lecture_13_enormous_exponents_and_card_shuffling.jl`)
    *   *Focus:* Efficient computation via successive squaring and modular invariants[cite: 1].
*   **Lecture 15: Open Secrets—Public Key Cryptography** (`lecture_15_open_secrets_public_key_cryptography.jl`)
    *   *Focus:* The RSA algorithm and the practical application of Euler's totient function[cite: 1].
*   **Lecture 16: The Birth of Graph Theory** (`lecture_16_the_birth_of_graph_theory.jl`)
    *   *Focus:* Core network topologies, Eulerian graphs, and Hamiltonian paths[cite: 1].
*   **Lecture 17: Ways to Walk—Matrices and Markov Chains** (`lecture_17_ways_to_walk_matrices_and_markov_chains.jl`)
    *   *Focus:* Transition probability matrices and stochastic modeling (e.g., PageRank)[cite: 1].
*   **Lecture 20: Weighted Graphs and Minimum Spanning Trees** (`lecture_20_weighted_graphs_and_minimum_spanning_trees.jl`)
    *   *Focus:* Greedy optimization algorithms for routing and connectivity[cite: 1].
*   **Lecture 23: Shortest Paths and Algorithm Complexity** (`lecture_23_shortest_paths_and_algorithm_complexity.jl`)
    *   *Focus:* Dijkstra’s Algorithm and computational complexity bounds (P vs. NP)[cite: 1, 2].

## The 3-Paradigm Methodology

Every mathematical theorem and proof within these notebooks is validated through a strict **Functional Programming & Data-Oriented Design (FP-DOD)** framework across three computational paradigms[cite: 1]:

1.  **Pure FP Focus (Math 1:1):** Pure, stateless data transformations using `map`, `reduce`, and `fold` without mutable state[cite: 1].
2.  **Abstract Algebra Focus:** Structural mathematics mapped into continuous algebraic geometry, polynomial rings, and Galois fields[cite: 1].
3.  **Group Theory & Geometric Focus:** Measurement of mathematical symmetries and topological spaces using the Symmetric Group and matrix projections[cite: 1].

## Technologies Used

*   **Julia:** The exclusive backend engine for all implementations, algorithms, and logic[cite: 1].
*   **Oscar.jl / AbstractAlgebra.jl:** Utilized for exact commutative algebra, ring bounds, and group theory operations[cite: 1].
*   **GraphRecipes.jl & Plots.jl:** Employed for visualizing graph topologies, minimum spanning trees, and network structures.
*   **Latexify.jl:** Used for parsing mathematical expressions and polynomial expansions into clean LaTeX formatting within outputs.
*   **Pluto.jl:** The reactive notebook environment housing all interactive lectures.