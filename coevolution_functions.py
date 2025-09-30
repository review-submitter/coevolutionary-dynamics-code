import numpy as np
from scipy.stats import truncnorm
import itertools
import csv


def get_truncnorm_number(mean, sd, low, up):
    """Get a random number that follow a truncated normal distribution.

    The truncated normal distribution is the probability distribution
    derived from that of a normally distributed random variable by
    bounding the random variable from either below or above (or both).

    Args:
        mean: The mean or expectation of the distribution.
        sd: The standard deviation of the distribution.
        low: The lower bound of the distribution.
        up: The upper bound of the distribution.

    Returns:
        A random number that follow a truncated normal distribution.
    """
    truncnorm_generator = truncnorm(
        (low - mean) / sd,
        (up - mean) / sd,
        loc=mean,
        scale=sd)
    truncnorm_random_number = truncnorm_generator.rvs()
    return truncnorm_random_number


def extinction_check(current_type, extinct_type, s=0):
    """Check extinction phenomenon.

    If one kind of prey or predator species go extinct, remove it from
    current_type to extinct_type list.

    Args:
        current_type: A list that stores the all the current type of prey
            or predator.
        extinct_type: A list that stores the extinct type of prey or predator.
        s: s=0. It used for counting the loop.
    """
    while s < len(current_type):
        if current_type[s].num == 0:
            extinct_type.append(current_type[s])
            del current_type[s]
        else:
            s += 1


def weighted_random(rng, a_weights):
    """Linear Scan"""
    weights = a_weights.flatten('C')
    remaining_distance = rng.random() * np.sum(weights)
    for i, weight in enumerate(weights):
        remaining_distance -= weight
        if remaining_distance < 0:
            if a_weights.ndim == 1:
                return i
            else:
                return i // a_weights.shape[1], i % a_weights.shape[1]


def find_reaction(rng, *reaction_rates):
    """ Select the reaction happen next time using Gillespie Algorithm
    and calculate the reaction time.

    Args:
        *reaction_rates: 1-D or 2-D array_like, array of reaction rate.

    Returns:
        a tuple of reaction time (tau), index of which series of
        reaction and index of reaction in that series.
    """
    # 1.Calculate total reaction rate:
    sum_reaction_rates_series = np.array([np.sum(i) for i in reaction_rates])
    total_reaction_rate = np.sum(sum_reaction_rates_series)

    # 2.Calculate the reaction time
    tau = - np.log(rng.random()) / total_reaction_rate

    # 3.Pick up the reaction to happen
    reaction_series_index = weighted_random(rng,sum_reaction_rates_series)
    reaction_index = weighted_random(rng, reaction_rates[reaction_series_index])

    return tau, reaction_series_index, reaction_index


def shannon_index(species_num_array):
    """Calculate the Shannon index: -∑pilnpi.

    The Shannon index is an information statistic index, which means it
    assumes all species are represented in a sample and that they are
    randomly sampled. p is the proportion (n/N) of individuals of one
    particular species found (n) divided by the total number of
    individuals found (N). The value of this index ranges between
    0 and 1, the greater the value, the greater the sample diversity.

    Args:
        species_num_array: An array that store the number of different kind
            of species.

    Returns:
        Shannon index of diversity of this population.
    """
    ratio = species_num_array / species_num_array.sum()
    shannon_index_diversity = - sum(ratio * np.log(ratio))
    return float('%0.4f' % shannon_index_diversity)


def simpson_index(species_num_array):
    """Calculate the Simpson's Diversity Index: 1 - ∑pi**2

    The Simpson index is a dominance index because it gives more weight to
    common or dominant species. In this case, a few rare species with only
    a few representatives will not affect the diversity.  p is the proportion
    (n/N) of individuals of one particular species found (n) divided by the
    total number of individuals found (N). The value of this index ranges
    between 0 and 1, the greater the value, the greater the sample
    diversity.

    Args:
        species_num_array: An array that store the number of different kind
            of species.

    Returns:
        Simpson's diversity index of this population.
    """
    ratio_ = species_num_array / species_num_array.sum()
    simpson_index_diversity = 1 - sum(ratio_ ** 2)
    return float('%0.4f' % simpson_index_diversity)


def list_append_data(list_list, data_list):
    i = 0
    while i < len(list_list):
        list_list[i].append(data_list[i])
        i += 1


def output_data(output_filename, written_data):
    with open(output_filename, 'a') as file_object:
        file_csv = csv.writer(file_object)
        file_csv.writerow(written_data)


# Function for getting truncated Laplace random variables:
def laplace_cdf(x, mu, b):
    """
    Compute the cumulative distribution function (CDF) of the Laplace distribution.

    Parameters:
    x (float or array-like): The input values at which to evaluate the CDF.
    mu (float): The mean (location) parameter of the Laplace distribution.
    b (float): The scale parameter (diversity of the distribution) of the Laplace distribution.

    Returns:
    float or array: The CDF values at the input points `x`.
    """
    # Standardize x by subtracting mu and dividing by scale parameter b
    z = (x - mu) / b
    # Apply the CDF formula: 0.5 * exp(z) for x < mu, 1 - 0.5 * exp(-z) for x >= mu
    return np.where(z < 0, 0.5 * np.exp(z), 1 - 0.5 * np.exp(-z))


def laplace_ppf(p, mu, b):
    """
    Compute the percent point function (PPF) or the inverse of the CDF for the Laplace distribution.

    Parameters:
    p (float or array-like): The probability values at which to evaluate the PPF.
    mu (float): The mean (location) parameter of the Laplace distribution.
    b (float): The scale parameter (diversity of the distribution) of the Laplace distribution.

    Returns:
    float or array: The corresponding quantiles (PPF values) for the given probabilities `p`.
    """
    # Convert p into a numpy array for element-wise operations
    p = np.asarray(p)
    
    # Conditional check: if p <= 0.5, use the formula for positive tail, otherwise for negative tail
    cond = p <= 0.5
    
    # Return quantiles based on the value of p
    return np.where(cond, mu + b * np.log(2 * p), mu - b * np.log(2 * (1 - p)))


def truncated_laplace_pdf(x, mu, b, x1, x2):
    """
    Compute the probability density function (PDF) of the truncated Laplace distribution.

    Parameters:
    x (float or array-like): The input values at which to evaluate the PDF.
    mu (float): The mean (location) parameter of the Laplace distribution.
    b (float): The scale parameter (diversity of the distribution) of the Laplace distribution.
    x1 (float): The lower bound of truncation.
    x2 (float): The upper bound of truncation.

    Returns:
    float or array: The PDF values at the input points `x` for the truncated Laplace distribution.
    """
    # Compute the CDF values at the truncation bounds x1 and x2
    F_x1 = laplace_cdf(x1, mu, b)
    F_x2 = laplace_cdf(x2, mu, b)
    
    # Normalize the PDF by dividing by the difference in CDF values (the normalization factor)
    normalization = F_x2 - F_x1
    
    # Create a mask for values of x within the truncation bounds [x1, x2]
    mask = (x >= x1) & (x <= x2)
    
    # Compute the PDF of the Laplace distribution and normalize it
    pdf = (1 / (2 * b)) * np.exp(-np.abs(x - mu) / b) / normalization
    
    # Apply the mask: return pdf if within bounds, otherwise 0
    return np.where(mask, pdf, 0.0)


def sample_truncated_laplace(mu, b, x1, x2):
    """
    Generate random samples from a truncated Laplace distribution.

    Parameters:
    mu (float): The mean (location) parameter of the Laplace distribution.
    b (float): The scale parameter (diversity of the distribution) of the Laplace distribution.
    x1 (float): The lower bound of truncation.
    x2 (float): The upper bound of truncation.

    Returns:
    array: An array of `n` samples drawn from the truncated Laplace distribution.
    """
    n=1
    # Compute the CDF values at the truncation bounds x1 and x2
    F_x1 = laplace_cdf(x1, mu, b)
    F_x2 = laplace_cdf(x2, mu, b)
    
    # Compute the difference in CDF values for normalization
    F_diff = F_x2 - F_x1
    
    # Generate uniform random numbers between 0 and 1
    u = np.random.uniform(size=n)
    
    # Transform the uniform random variables to lie within the truncated CDF range
    u_truncated = F_x1 + u * F_diff
    
    # Use the inverse CDF (PPF) to obtain the corresponding Laplace samples
    samples = laplace_ppf(u_truncated, mu, b)
    
    return samples[0]

