from utils import get_distance_matrix
from utils import get_solution_cost
from utils import get_swap_neighbourhood

class ProblemInstance:
    def __init__(self, n, p, alpha, delta, ksi, node_coordinates, demand):
        self.n = n
        self.p = p
        self.alpha = alpha
        self.delta = delta
        self.ksi = ksi
        self.node_coordinates = node_coordinates
        self.demand = demand
        self.distances = get_distance_matrix(self.node_coordinates)

class Solution:
    def __init__(self, hubs, problem, cost=None):
        self.hubs = hubs
        self.problem = problem
        if cost is not None:
            self.cost = cost
        else:
            self.cost = self._get_cost()
        # list of hubs lists (not list of solutions)
        self.swap_neighbourhood = get_swap_neighbourhood(hubs, problem.n)

    def __str__(self):
        return f"Solution(hubs={self.hubs})"

    def _get_cost(self):
        return get_solution_cost(self.hubs, self.problem)

    def get_neighbourhood(self, neighbourhood_type):
        if neighbourhood_type == 'swap':
            return self.swap_neighbourhood
        else:
            raise ValueError(f"Unknown neighbourhood type. Supported neighbourhood types: {NEIGHBOURHOOD_TYPES}")
