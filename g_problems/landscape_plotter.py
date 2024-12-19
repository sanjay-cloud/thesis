import numpy as np
import matplotlib.pyplot as plt


class plotter:
    def __init__(self, problem, fig_length = 5, fig_width = 4, x_ticks = 1, discrete_steps = 1, x_ticks_rotate = False):
        self.problem = problem
        self._integer_var = self.integer_var()
        self._real_var = self.real_var()
        self.fig_length = fig_length
        self.fig_width = fig_width
        self.x_ticks = x_ticks
        self.x_ticks_rotate = x_ticks_rotate 
        self.discrete_steps = discrete_steps

    def mix_int_prob(self,discrete_val, continuous_val):

        xopt_copy = self.problem.xopt().copy()
        xopt_copy = np.array(xopt_copy)
        integer_index = self.problem.integer_indices()[0]
        real_index = self.problem.real_indices()[0]
        integer_indices = self.problem.integer_indices()
        xopt_copy[integer_index] = discrete_val
        xopt_copy[real_index] = continuous_val
        xopt_copy[integer_indices] = np.round(xopt_copy[integer_indices])

        return self.problem.objective(xopt_copy)
    
    def integer_var(self, index = 0):
        return self.problem.integer_indices()[index]
    
    def real_var(self, index = 0):
        return self.problem.real_indices()[index]
    
    def real_bounds(self):
        return self.problem.lower_bounds()[self._real_var], self.problem.upper_bounds()[self._real_var]
    
    def integer_bounds(self):
        return self.problem.lower_bounds()[self._integer_var], self.problem.upper_bounds()[self._integer_var]
        
    
    def get_results(self, discrete_addition = 0, discrete_steps = 1, mixed_integer = True):
        dis_lower_bounds, dis_upper_bounds = self.integer_bounds()
        con_lower_bounds, con_upper_bounds = self.real_bounds()

        # step_size = get_discrete_step_size(dis_lower_bounds, dis_upper_bounds)
        if mixed_integer:
            continuous_var_range = np.linspace(con_lower_bounds, con_upper_bounds, 100) 
        else:
            continuous_var_range = np.arange(con_lower_bounds + discrete_addition, con_upper_bounds + 1, discrete_steps)    
        discrete_var_range = np.arange(dis_lower_bounds + discrete_addition, dis_upper_bounds + 1, discrete_steps) 
        results = np.zeros((len(discrete_var_range), len(continuous_var_range)))
        
        for i, discrete_val in enumerate(discrete_var_range):  
            if mixed_integer:          
                continuous_var_range = np.linspace(con_lower_bounds, con_upper_bounds, 100)  
            else:   
                continuous_var_range = np.arange(con_lower_bounds + discrete_addition , con_upper_bounds + 1, discrete_steps) 
            for j, continuous_val in enumerate(continuous_var_range):
                results[i, j] = self.mix_int_prob(discrete_val, continuous_val)
        return  np.meshgrid(continuous_var_range, discrete_var_range) , results
    
    def plot_heatmap(self):
        if self.problem.problem_type:
            y_label = "Continuous"
        else :
            y_label = "Discrete"  
        dis_lower_bounds, dis_upper_bounds = self.integer_bounds()
        (X0, X1), Z = self.get_results(mixed_integer=self.problem.problem_type)  
        plt.figure(figsize=(self.fig_length, self.fig_width))
        plt.pcolormesh(X1, X0, Z, shading='auto', cmap='viridis')
        plt.colorbar(label='Objective Value')
        plt.xlabel('x'+str(self._integer_var+1)+'(Discrete)')
        plt.ylabel('x'+str(self._real_var+1)+' ('+y_label+')')
        plt.title('Heatmap of'+self.problem.fn())
        plt.xticks(np.arange(dis_lower_bounds , dis_upper_bounds + 1, self.x_ticks) )
        if self.x_ticks_rotate:
            plt.xticks(rotation=90)


    def plot_level_sets(self):
        
        if self.problem.problem_type:
            y_label = "Continuous"
        else :
            y_label = "Discrete" 
        (X, Y), results = self.get_results(discrete_addition = 0.5, mixed_integer=self.discrete_steps)
        plt.figure(figsize=(self.fig_length, self.fig_width))
        cp = plt.contour(Y, X, results, levels=20, cmap='viridis')

        plt.clabel(cp, inline=True, fontsize=10, fmt='%1.1f')  # Add labels to the contour lines
        plt.colorbar(cp)
        plt.xlabel('Discrete variable (x'+str(self._integer_var+1)+')')
        plt.ylabel(y_label+' variable (x'+str(self._real_var+1)+')')        
        plt.title(self.problem.fn())
        plt.show()