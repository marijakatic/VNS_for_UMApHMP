from utils import get_nodes
from utils import Solution

from operator import attrgetter
import random
from tqdm import tqdm

from global_parameters import NEIGHBOURHOOD_TYPES


MAX_ITER = 25
PRECISION = 0.0001

def get_initial_solution(n, p, distances):
    """Solution initialization method for vns.

    Paramteres: 
    n (int): Number of nodes
    p (int): Number of hubs
    distances (matrix): distance matrix for the graph

    Returns:
    (list): hubs

    Notes:
    This initialization method is taken from the paper
    El-Ghazali Talbi, Raca Todosijevic, 'The robust uncapacitated multiple allocation p-hub median problem' 

    "A solution built in a greedy manner. 
    Namely, the initial p hubs are chosen as those whose maximum transportation cost
    to any other node are the p smallest. More precisely, 
    let g(p) be the maximum transportation cost of a node h to any other node, 
    i.e. g(h) = max{C_ih|i from N, i != h}, h from N.
    Then, the nodes with the p smallest values of function g are taken as the initial p hubs."
    """
    nodes = get_nodes(n)
    longest_edge_from_node = {}
    for node in nodes:
        longest_edge_from_node[node] = max(distances[node])
    hubs = sorted(longest_edge_from_node, key=longest_edge_from_node.get)[:p]
    return hubs

def get_best_solution(solutions):
    return min(solutions, key=attrgetter('cost'))

# Lets use Best Improvement Search. We can also try First Improvement search later.
def local_search(solution, neighbourhood_type):
    ''' Best Improvement Search'''
    neighbourhood = [Solution(neighbour, solution.problem)
                     for neighbour in solution.get_neighbourhood(neighbourhood_type)]
    while True:
        curr_solution = solution
        # solution ← argmin{f(s)}, s ∈ N(solution)
        solution = get_best_solution(neighbourhood)
        # if no direction of descent anymore
        if solution.cost >= curr_solution.cost:
            break
    return curr_solution

def shake(solution, neighbourhood_type):
    # print(f"shaking solution={solution}, neighbourhood={solution.get_neighbourhood(neighbourhood_type)}")
    return Solution(random.choice(solution.get_neighbourhood(neighbourhood_type)), 
                    solution.problem)


def basic_VNS(problem, max_iter=MAX_ITER, precision=PRECISION, verbose=False):
    # initialize solution
    solution = Solution(get_initial_solution(problem.n, problem.p, problem.distances), problem)
    optimal_solution = solution
    if verbose == True:
        iter_range = tqdm(range(max_iter))
    else: 
        iter_range = range(max_iter)
    for iteration in iter_range:
        i = 0
        while i < len(NEIGHBOURHOOD_TYPES):
            # Shaking
            rand_solution = shake(solution, NEIGHBOURHOOD_TYPES[i])
            # print("shake")
            # print(f"rand_solution={rand_solution}")
            # Local search
            local_min = local_search(rand_solution, NEIGHBOURHOOD_TYPES[i])
            # print("local_search")
            # print(f"local_min={local_min}")
            # Change neighbourhood
            if local_min.cost < solution.cost:
                solution = local_min
                # reset loop 
                i = 0
            else:
                i += 1
            # print(i)
        # Update optimal solution
        if solution.cost < optimal_solution.cost:
            optimal_solution = solution
        if problem.optimal_cost is not None \
            and abs(optimal_solution.cost - problem.optimal_cost) < PRECISION:
            break
    return optimal_solution