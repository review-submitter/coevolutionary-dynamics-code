class Prey:
    """
    A class to represent a type of prey species.

    Attributes:
        num: Number of this prey species.
        g_value: Growth-defense trade off parameter, g ∈ (0,1).
        num_list: A list that store the number of this species after
        each reaction happen.
        time_list: A list that store the time(One-to-one correspondence
        to the num_list).
    """

    def __init__(self, num, g_value):
        """Initialize attributes to describe a kind of prey."""
        self.num = num
        self.g_value = float('%.3f' % g_value)
        self.num_list = []
        self.time_list = []

    def prey_reproduce(self):
        """Prey reproduces without mutation."""
        self.num += 1

    def prey_die(self):
        """Prey dies because of competition, predation or dies
        naturally."""
        self.num -= 1

    @staticmethod
    def prey_mutation(random_distribution_function,
                      current_prey_list, *truncnorm_args):
        """Mutation of parental prey to generate a new prey species.

        When a prey reproduces, a mutation may occur with a small
        probability μx. The mutant is characterized by a new g value
        (type specific rate) drawn randomly from uniform distribution
        or truncated normal distribution between 0 and 1. (In truncated
        normal distribution, the mean of this distribution is g value
        of the parental prey. The lower bound is 0 and the upper bound
        is 1.)

        Args:
            random_distribution_function: "uniform" or "Gaussian", determine
            which random distribution that the new g value(type specific
            rate) drawn from.
            current_prey_list: a list that store the current prey type.
            *truncnorm_args: the standard deviation and other parameters
            of Gaussian distribution. If the mutant prey's g value is drawn
            from truncated normal distribution, this argument is needed.
        """
        new_prey = Prey(num=1,
                        g_value=random_distribution_function(*truncnorm_args))
        current_prey_list.append(new_prey)

    def record_data(self, time_):
        """Record the number of prey and time after each reaction happen."""
        self.time_list.append(time_)
        self.num_list.append(self.num)


class Predator:
    """
    A class to represent a type of predator species.

    Attributes:
        num: Number of this predator species.
        k_value: k value is the ratio of predator growth to predation
            and represents the reproduction efficacy of a predator type.
        num_list: A list that store the number of this species after a
            reaction happen.
        time_list: A list that store the time.
    """

    def __init__(self, num, k_value):
        """Initialize attributes to describe a kind of predator."""
        self.k_value = float('%.3f' % k_value)
        self.num = num
        self.num_list = []
        self.time_list = []

    def predator_reproduce(self):
        """Predator reproduces without mutation."""
        self.num += 1

    def predator_die(self):
        """Predator dies naturally."""
        self.num -= 1

    @staticmethod
    def predator_mutation(random_distribution_function,
                          current_predator_list, *truncnorm_args):
        """Mutation of parental predator to generate a new prey species.

        With a probability μy, a predator produces a mutant with a new
        k value drawn from a uniform distribution or truncated normal
        distribution between 0 and kmax.
        (In truncated normal distribution:the upper limit of the
        reproduction efficacy. The mean of this distribution is the k
        value of the parental predator. The lower bound is 0 and the
        upper bound is 0.3.

        Args:
            random_distribution_function: "uniform" or "Gaussian", determine
            which random distribution that the new k value drawn from.
            current_predator_list: a list that store the current
                predator type.
            *truncnorm_args: the standard deviation and other parameters
            of Gaussian distribution, if the mutant predator's k value is
            drawn from truncated normal distribution, this argument is needed.
        """
        new_predator = Predator(num=1,
                                k_value=random_distribution_function(
                                    *truncnorm_args))
        current_predator_list.append(new_predator)

    def record_data(self, time_):
        """Record the number of prey and time after each reaction happen."""
        self.time_list.append(time_)
        self.num_list.append(self.num)
