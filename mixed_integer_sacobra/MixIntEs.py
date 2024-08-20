import numpy as np

class MixedIntegerOptimizer:
    def __init__(self, obj_func, max_eval=np.inf, lower=None, upper=None, minimize=True, mu=50, lambda_=10, integer_indices = [], verbose=False, constraints_fun = None):
        self.obj_func = obj_func  # Objective function
        self.max_eval = max_eval  # Maximum evaluations
        self.lower = np.array(lower)  # Lower bounds
        self.upper = np.array(upper)  # Upper bounds
        self.mu = mu  # Population size
        self.lambda_ = lambda_  # Number of offspring
        self.eval_count = 0  # Counter for evaluations
        self.iter_count = 0  # Counter for iterations
        self.minimize = minimize  # Optimization goal
        self.verbose = verbose  # Verbosity flag
        self.curr_gen = 1  # Current generation
        self.integer_indices = integer_indices # Indices of integer variables
        # self.initial_pop = self.generate_pop(pop_size=self.mu)
        self.constraints_fun = constraints_fun

    def generate_pop(self, pop_size, integer_indices):
        pop = np.random.uniform(low=self.lower, high=self.upper, size=(pop_size, len(self.lower)))
        pop[:, self.integer_indices] = np.round(pop[:, self.integer_indices])
        return pop

    def recombine(self, parent1, parent2):
        child = np.copy(parent1)

        crossover_mask = np.random.rand(len(child)) > 0.5
        child = np.where(crossover_mask, parent1, parent2)
        return child

    def mutate(self, individual):
        # Mutate continuous variables
        for i in range(len(individual)):
            if i not in self.integer_indices:
                individual[i] += np.random.normal(0, 0.1)  # Gaussian mutation

        # Mutate discrete variables
        for i in self.integer_indices:
            if np.random.rand() < 0.5:  # Mutation probability
                individual[i] = np.round(np.random.uniform(self.lower[i], self.upper[i]))

        # Enforce bounds
        individual = np.clip(individual, self.lower, self.upper)
        return individual

    def dynamic_penalty(self, X, t, C=0.5, alpha=1, minimize=True):
        p = (-1) ** (not minimize) * (C * 1) ** alpha * np.sum(X)
        return p

    def evaluate(self, population):
        fitness = np.zeros(len(population))
        for i, individual in enumerate(population):
            # obj_value, *constraints = self.obj_func(individual)
            # penalty = self.dynamic_penalty(constraints, self.curr_gen)
            obj_value = self.obj_func(individual)
            constraints_val = self.constraints_fun(individual)
            penalty = self.dynamic_penalty(constraints_val, self.curr_gen)
            fitness[i] = obj_value + penalty if self.minimize else -obj_value + penalty
            # obj_value = self.obj_func(individual)
            # fitness[i] = -obj_value  if self.minimize else -obj_value 
        self.eval_count += len(population)
        return fitness

    def select(self, pop, fitness):
        rank = np.argsort(fitness)
        return pop[rank[:self.mu]], fitness[rank[:self.mu]]

    def stop(self):
        return self.eval_count >= self.max_eval

    def optimize(self, pop):
        # pop = self.initial_pop
        fitness = self.evaluate(pop)

        while not self.stop():
            offspring = []
            for _ in range(self.lambda_):
                parents = np.random.choice(range(self.mu), 2, replace=False)
                child = self.recombine(pop[parents[0]], pop[parents[1]])
                child = self.mutate(child)
                offspring.append(child)

            offspring = np.array(offspring)
            fitness_offspring = self.evaluate(offspring)
            pop, fitness = self.select(np.vstack((pop, offspring)), np.concatenate((fitness, fitness_offspring)))

            if self.verbose:
                best_fitness = fitness[0]
                print(f'Iteration {self.iter_count + 1}, Best Fitness: {best_fitness}, Best Individual: {pop[0]}')

            self.iter_count += 1
            self.curr_gen += 1

        return pop[0], fitness[0]
def G01(x):
    obj = np.sum(5*x[:4])-(5*np.sum(x[:4]**2))-(np.sum(x[4:13]))
    g1 = 2*x[0]+2*x[1]+x[9]+x[10] - 10
    g2 = 2*x[0]+2*x[2]+x[9]+x[11] - 10
    g3 = 2*x[1]+2*x[2]+x[10]+x[11] - 10
    
    g4 = -8*x[0]+x[9]
    g5 = -8*x[1]+x[10]
    g6 = -8*x[2]+x[11]
    
    g7 = -2*x[3]-x[4]+x[9]
    g8 = -2*x[5]-x[6]+x[10]
    g9 = -2*x[7]-x[8]+x[11]
    
    return(np.array([obj,g1,g2,g3,g4,g5,g6,g7,g8,g9]))

def G05(x):
    obj = 3 * x[0] + 0.000001 * x[0]**3 + 2 * x[1] + (0.000002 / 3) * x[1]**3
    # return obj
    g1 = -x[3] + x[2] - 0.55
    g2 = -x[2] + x[3] - 0.55

    h31 = -1 * (1000 * np.sin(-x[2] - 0.25) + 1000 * np.sin(-x[3] - 0.25) + 894.8 - x[0] + 0.0001)
    h32 = 1000 * np.sin(-x[2] - 0.25) + 1000 * np.sin(-x[3] - 0.25) + 894.8 - x[0] - 0.0001
    h41 = -1 * (1000 * np.sin(x[2] - 0.25) + 1000 * np.sin(x[2] - x[3] - 0.25) + 894.8 - x[1] + 0.0001)
    h42 = 1000 * np.sin(x[2] - 0.25) + 1000 * np.sin(x[2] - x[3] - 0.25) + 894.8 - x[1] - 0.0001
    h51 = -1 * (1000 * np.sin(x[3] - 0.25) + 1000 * np.sin(x[3] - x[2] - 0.25) + 1294.8 + 0.0001)
    h52 = 1000 * np.sin(x[3] - 0.25) + 1000 * np.sin(x[3] - x[2] - 0.25) + 1294.8 - 0.0001

    return np.array([obj, g1, g2, h31, h32, h41, h42, h51, h52])
def run_05():
    lower_bounds = np.array([0, 0, -0.55, -0.55])
    upper_bounds = np.array([1200, 1200, 0.55, 0.55])
    integer_indices = [1]
    optimizer = MixedIntegerOptimizer(obj_func=G05, lower=lower_bounds, upper=upper_bounds, max_eval=1000, minimize=True, verbose=True, integer_indices = integer_indices)

    initial_population = optimizer.generate_pop(pop_size=10, integer_indices=integer_indices)
    print("initial pop", initial_population[0])
    best_solution, best_fitness = optimizer.optimize(initial_population)

    print(f'Best Solution: {best_solution}, Best Fitness: {best_fitness}')

def G02(x):
    obj = (1-x[0])**2 + 100*(x[1]-x[0]**2)**2
    g1 = (x[0]-1)**3-x[1]+1
    g2 = x[0]+x[1]-2
    res = np.array([obj,g1,g2])
    return(obj)


def run_02():
    lower_bounds = np.array([-1.5, -0.5])
    upper_bounds = np.array([1.5, 2.5])
    # nConstraints=2
    # d=4
    integer_indices = [1]
    optimizer = MixedIntegerOptimizer(obj_func=G02, lower=lower_bounds, upper=upper_bounds, max_eval=1000, minimize=True, verbose=True, integer_indices = integer_indices)

    initial_population = optimizer.generate_pop(pop_size=10, integer_indices=integer_indices)
    print("initial pop", initial_population[0])
    best_solution, best_fitness = optimizer.optimize(initial_population)

    print(f'Best Solution: {best_solution}, Best Fitness: {best_fitness}')

# run_05()    