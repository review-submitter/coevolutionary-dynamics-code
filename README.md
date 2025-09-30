# Predator-Prey Coevolution Simulation

This repository contains code for simulating predator-prey coevolution dynamics using stochastic individual-based models with the Gillespie algorithm. The simulation models the evolutionary dynamics of prey and predator populations with different mutation distribution models. 

## Files

- **`species.py`**: Defines `Prey` and `Predator` classes with reproduction, death, and mutation methods
- **`coevolution_functions.py`**: Utility functions including Gillespie algorithm, distribution samplers, and extinction checks
- **`model_truncated_normal_distribution.py`**: Simulation with truncated normal mutation distribution
- **`model_truncated_laplace_distribution.py`**: Simulation with truncated Laplace mutation distribution
- **`model_uniform_distribution.py`**: Simulation with uniform mutation distribution

## Required Software and Packages

- Python 3.7 or higher  
- numpy  
- scipy  
- pandas
- 
## Usage

### Truncated Normal/Laplace Models

```bash
python model_truncated_normal_distribution.py <m> <std_dev> <ancestor_prey_g> <ancestor_pred_k> <rep> <seed>
```

### Uniform Model

```bash
python model_uniform_distribution.py <m> <ancestor_prey_g> <ancestor_pred_k> <rep> <seed>
```

**Arguments:**

- `m`: Predation scaling coefficient (e.g., 0.5, 1.0, 2.0)
- `std_dev`: Mutation standard deviation (e.g., 0.1, 0.2)
- `ancestor_prey_g`: Initial prey trait value, range [0,1] (e.g., 0.5)
- `ancestor_pred_k`: Initial predator trait value, range [0,0.3] (e.g., 0.15)
- `rep`: Replication ID (e.g., 1, 2, 3)
- `seed`: Random seed (e.g., 12345)

**Example:**

```bash
python model_truncated_normal_distribution.py 1.0 0.1 1.0 0.3 1 12345
```

## Output

Generates `.pkl` files containing:

- `status`: 0 (both extinct), 1 (predator extinct), or 2 (coexist)
- `simulation_time`: Time when simulation ended
- `output_dict`: DataFrame with time series of trait values and population sizes

**Read output:**

```python
import pickle
with open('simulation_data_1.pkl', 'rb') as f:
    data = pickle.load(f)
```

## Key Parameters

- Initial populations: 1000 prey, 100 predators
- Max simulation time: 2500 units
- Prey mutation rate: 0.0005
- Predator mutation rate: 0.005
- Recording interval: Every 1 time unit

## Note

File paths are parameterized for anonymity. To customize output location, modify the `output_path` variable in the script's main block.
