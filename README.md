# Predator-Prey Coevolution Simulation
This repository contains code for simulating predator-prey coevolution dynamics using stochastic individual-based models with the Gillespie algorithm. The simulation models the evolutionary dynamics of prey and predator populations with different mutation distribution models. 

## Overview
The simulation tracks the coevolution of prey and predator species where:
- **Prey** have a growth-defense trade-off parameter (g-value, ranging from 0 to 1)
- **Predator** have a reproduction efficiency parameter (k-value, ranging from 0 to 0.3)
- Mutations occur during reproduction, introducing trait variation
- Population dynamics are governed by birth, death, predation, and competition processes


Files Description

### Core Files

1. **`species.py`**
    - Defines the `Prey` and `Predator` classes
    - Implements methods for reproduction, death, and mutation
    - Each species type is characterized by its trait value and population size
2. **`coevolution_functions.py`**
    - Contains utility functions for the simulation:
        - `get_truncnorm_number()`: Generates random numbers from truncated normal distribution
        - `sample_truncated_laplace()`: Generates random numbers from truncated Laplace distribution
        - `find_reaction()`: Implements Gillespie algorithm for stochastic simulation
        - `extinction_check()`: Checks and handles species extinction
        - Diversity indices: `shannon_index()` and `simpson_index()`
        - Additional helper functions for data recording
3. **`model_truncated_normal_distribution.py`**
    - Main simulation script using **truncated normal distribution** for mutations
    - Mutations are drawn from a normal distribution centered at parent trait value
4. **`model_truncated_laplace_distribution.py`**
    - Main simulation script using **truncated Laplace distribution** for mutations
    - Laplace distribution has heavier tails than normal distribution
5. **`model_uniform_distribution.py`**
    - Main simulation script using **uniform distribution** for mutations
    - Mutations are uniformly distributed across the entire trait range
