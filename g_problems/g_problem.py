from landscape_plotter import plotter as plt
import numpy as np
class GProblem():
    def __init__(self, fn, dimensions, objective, constraints, lower_bounds, upper_bounds, integer_indices, real_indices, minimize, xopt, fopt, problem_type_mixint = True, x_ticks = 1, x_ticks_rotate = False, discrete_steps = 1,  plotter = None):
        self._fn = fn
        self._dimensions = dimensions
        self.problem_type = problem_type_mixint
        self._objective = objective
        self._constraints = constraints
        # self._combined_obj = self.combined_obj
        self._lower_bounds = lower_bounds
        self._upper_bounds = upper_bounds
        self._real_indices = real_indices
        self._integer_indices = integer_indices
        self.minimize = minimize
        self._fopt = fopt
        self._xopt = xopt
        # self.x_ticks = x_ticks
        self.plotter = plt(self, x_ticks = x_ticks, discrete_steps = discrete_steps, x_ticks_rotate = x_ticks_rotate) if plotter == None else plotter

    def fn(self):
        return str(self._fn)

    def dimensions(self):
        return self._dimensions   

    def constraints(self, x):
        # x_discrete = self.discretize_variables(x)
        return self._constraints(x)

    def objective(self, x):
        # x_discrete = self.discretize_variables(x)
        return self._objective(x)  if self.minimize  else -self._objective(x)

    def lower_bounds(self):
        return self._lower_bounds

    def upper_bounds(self):
        return self._upper_bounds   
    
    def combined_obj(self, x):
        return np.array([self.objective(x)] + (self.constraints(x)))
    
    def real_indices(self):
        return self._real_indices
    
    def discretized_x(self, x):
        return self.discretize_variables(x)
    
    def integer_indices(self):
        return self._integer_indices
    
    def xopt(self):
        return self._xopt
    
    def fopt(self):
        return self._fopt
    
    # Example: Plotting the G01 problem
    def plot_landscape(self, heatmap = True):
        if heatmap:
            self.plotter.plot_heatmap()
        else:    
            self.plotter.plot_level_sets()
