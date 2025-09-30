import logging
import argparse
import numpy as np
import pandas as pd
from species import Prey, Predator
import coevolution_functions as cf
import os
import pickle

# Setting up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class Parameters:
    """
    Class to hold the global parameters for the model.
    """
    def __init__(self):
        self.prey_growth_rate = 1.0  # baseline growth rate of prey
        self.prey_death_rate = 0.1  # intrinsic death rate of prey
        self.predator_death_rate = 0.5  # intrinsic death rate of predator
        self.resource_competition_coef = 0.00005  # resource competition coefficient
        self.prey_mutation_rate = 0.0005  # mutation rate of prey per division
        self.predator_mutation_rate = 0.005  # mutation rate of predation per division
        self.predation_scaling_coef = 0.005  # scaling coefficient of the predation rate
        self.max_predation_capacity = 0.3  # maximum predation capacity


def initialize_parameters():
    """
    Initialize and return the global parameters for the model.
    """
    return Parameters()


def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description='Run coevolution simulation.')
    parser.add_argument('m', type=float, help='m determines the initial shape of the growth-defense trade-off')
    parser.add_argument('ancestor_prey_growth', type=float, help='Ancestor prey growth value')
    parser.add_argument('ancestor_predator_capacity', type=float, help='Ancestor predator capacity value')
    parser.add_argument('rep', type=int, help='Array job replication ID')
    parser.add_argument('seed', type=int, help='Random number generator seed')
    parser.add_argument('--output_dir', type=str, default='./output', 
                   help='Output directory path')

    args = parser.parse_args()
    return (args.m, args.ancestor_prey_growth,
            args.ancestor_predator_capacity, args.rep, args.seed)


def record_simulation_data(output_filename, status, simulation_time, output_dict):
    """
    Record the simulation data to the output file using pickle.
    """
    data = {
        'status': status,
        'simulation_time': simulation_time,
        'output_dict': output_dict
    }
    
    with open(output_filename, 'wb') as f:
        pickle.dump(data, f)

def record_dict(output_dict, time_index,time, prey_type, prey_pop, predator_type, predator_pop):
    new_row = {'time_index': time_index,
               'time': time,
               'prey_type': prey_type,
               'prey_population': prey_pop,
               'predator_type': predator_type,
                'predator_population': predator_pop}
    new_row_df = pd.DataFrame([new_row])
    output_dict = pd.concat([output_dict, new_row_df], ignore_index=True)
    return output_dict

def run_simulation(params,m, ancestor_prey_growth,
    ancestor_predator_capacity, output_filename, seed):
    """
    Run the coevolution simulation.
    """  
    
    rng = np.random.default_rng(seed)
    
    logging_info = (
        f"Simulation Begins! Uniform model! Predation scaling coefficient m={m}, "
        f"Ancestor prey growth={ancestor_prey_growth}, "
        f"Ancestor predator capacity={ancestor_predator_capacity}"
    )
    logging.info(logging_info)
    
    # set up initial conditions

    ancestor_prey = Prey(num=1000, g_value=ancestor_prey_growth)
    ancestor_predator = Predator(num=100, k_value=ancestor_predator_capacity)
    
    t = np.float64(0)
    max_time = np.float64(2500)


    current_prey = [ancestor_prey]
    current_predator = [ancestor_predator]
    extinct_prey = []
    extinct_predator = []
    
    output_dict = pd.DataFrame(columns=['time_index', 'time','prey_type', 'prey_population', 'predator_type', 'predator_population'])
    
    time_list = np.arange(0, 2501, 1)
    time_index = 0


    while t <= max_time:
        cf.extinction_check(current_prey, extinct_prey)
        cf.extinction_check(current_predator, extinct_predator)
        
        prey_growth_values = np.array([prey.g_value for prey in current_prey])
        prey_population = np.array([prey.num for prey in current_prey])
        predator_capacity_values = np.array([predator.k_value for predator in current_predator])
        predator_population = np.array([predator.num for predator in current_predator])
        
        if t >= time_list[time_index]:
            output_dict = record_dict(output_dict, time_index,t,
                                      prey_growth_values,
                                      prey_population,
                                      predator_capacity_values,
                                      predator_population)
            
            time_index +=1
    
        if len(current_prey) == 0:
            output_dict = record_dict(output_dict, time_index,t,
                                      prey_growth_values,
                                      prey_population,
                                      predator_capacity_values,
                                      predator_population)
            
            record_simulation_data(output_filename, 0, t,
                output_dict)
                
            logging.info(f"Time {t}: All prey extinct. Both species extinct!")
            break

        elif len(current_predator) == 0:
            output_dict = record_dict(output_dict, time_index,t,
                                      prey_growth_values,
                                      prey_population,
                                      predator_capacity_values,
                                      predator_population)
            
            record_simulation_data(output_filename, 1, t,output_dict
               )

            logging.info(f"Time {t}: Predators extinct!")
            break

        else:
            prey_birth_no_mutation = params.prey_growth_rate * (1 - params.prey_mutation_rate) * prey_population * prey_growth_values
            prey_birth_with_mutation = params.prey_growth_rate * params.prey_mutation_rate * prey_growth_values * prey_population
            prey_competition_death = (prey_population[:, np.newaxis] * prey_population) * params.resource_competition_coef
            prey_intrinsic_death = params.prey_death_rate * prey_population

            predation_rate = params.predation_scaling_coef * prey_growth_values[:, np.newaxis] ** (m*(
                predator_capacity_values / params.max_predation_capacity))

            predation_no_birth = predator_population * (1 - predator_capacity_values) * prey_population[:, np.newaxis] * predation_rate
            predation_birth_no_mutation = predator_population * (1 - params.predator_mutation_rate) * predator_capacity_values * prey_population[:, np.newaxis] * predation_rate
            predation_birth_with_mutation = predator_population * predator_capacity_values * params.predator_mutation_rate * prey_population[:, np.newaxis] * predation_rate
            predator_intrinsic_death = params.predator_death_rate * predator_population


            tau, reaction_series_index, reaction_index = cf.find_reaction(
                rng,
                prey_birth_no_mutation, prey_birth_with_mutation,
                prey_competition_death, prey_intrinsic_death,
                predation_no_birth, predation_birth_no_mutation,
                predation_birth_with_mutation, predator_intrinsic_death)

            t += np.float64(tau)

            if reaction_series_index == 0:
                current_prey[reaction_index].prey_reproduce()
            elif reaction_series_index == 1:
                current_prey[reaction_index].prey_mutation(rng.uniform,current_prey, 0, 1.0)
            elif reaction_series_index == 2:
                current_prey[reaction_index[0]].prey_die()
            elif reaction_series_index == 3:
                current_prey[reaction_index].prey_die()
            elif reaction_series_index == 4:
                current_prey[reaction_index[0]].prey_die()
            elif reaction_series_index == 5:
                current_prey[reaction_index[0]].prey_die()
                current_predator[reaction_index[1]].predator_reproduce()
            elif reaction_series_index == 6:
                current_prey[reaction_index[0]].prey_die()
                current_predator[reaction_index[1]].predator_mutation(rng.uniform, current_predator, 0, 0.3)
            elif reaction_series_index == 7:
                current_predator[reaction_index].predator_die()


    if len(current_prey) != 0 and len(current_predator) != 0:
        logging.info(f"Time {t}: Both species coexist!")
        
        # save the output_prey_list, output_predator_list, output_prey_population_list, output_predator_population_list to the output_filename
        record_simulation_data(output_filename, 2, t, output_dict)

    logging.info('Simulation end!')

if __name__ == "__main__":
    params = initialize_parameters()
    args = parse_arguments()
    
    output_path = os.path.join(args.output_dir, 'uniform', 
                          f'm={args[0]}}')
    # if the output path does not exist, create it
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)
    output_filename = output_path + 'simulation_data_{}.pkl'.format(args[3])
    run_simulation(params, *args[:3], output_filename, args[4])
