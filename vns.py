from utils import get_nodes
from utils import Solution

from operator import attrgetter
import random
from tqdm import tqdm
import time
import math
import numpy as np

from utils import NEIGHBOURHOOD_TYPES

T_MAX = 2000 # in seconds
MAX_ITER = 20
PRECISION = 0.01

def get_initial_solution_random(n, p, _, nodes_coor=None):
    return random.sample(range(n), p)

def get_initial_solution_circle(n, p, _, nodes_coor):
    center_x, center_y = np.average([x for x,_ in nodes_coor]),  np.average([y for _,y in nodes_coor])
    nodes_coor = [(x - center_x, y - center_y) for x,y in nodes_coor]

    angle = (2*math.pi)/p
    bucket_angles = {hub:(hub*angle, (hub+1)*angle) for hub in range(p)}
    buckets = {hub:[] for hub in range(p)}

    for i in range(n):
        x,y = nodes_coor[i]
        node_angle = math.atan(y/x)
        if x<0 and y>0:
            node_angle += math.pi
        if x<0 and y<0:
            node_angle += math.pi
        for bucket in bucket_angles:
            if node_angle >= bucket_angles[bucket][0] and node_angle < bucket_angles[bucket][1]:
                buckets[bucket].append(i)

    hubs = []
    for bucket in buckets:
        if len(buckets[bucket]) > 0:
            hubs.append(random.choice(buckets[bucket]))

    hubs.extend(random.sample([i for i in range(n) if i not in hubs], p-len(hubs)))
    return hubs

def get_initial_solution_robust(n, p, distances, nodes_coor=None):
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
    # solution ← argmin{f(s)}, s ∈ solutions
    return min(solutions, key=attrgetter('cost'))

def local_search_best_improvement(solution, neighbourhood_type, use_c):
    ''' Best Improvement Search'''
    while True:
        curr_solution = solution
        neighbourhood = [Solution(neighbour, solution.problem, use_c=use_c)
                     for neighbour in solution.get_neighbourhood(neighbourhood_type)]
        solution = get_best_solution(neighbourhood)
        # if no direction of descent anymore
        if solution.cost >= curr_solution.cost:
            break
    return curr_solution

def local_search_first_improvement(solution, neighbourhood_type, use_c):
    '''First Improvement Search'''
    while True:
        curr_solution = solution
        # take the first better solution, not the best, from neighbourhood
        neighbourhood = solution.get_neighbourhood(neighbourhood_type)
        for neighbour in neighbourhood:
            neighbour = Solution(neighbour, solution.problem, use_c)
            # early stopping, when first descent is found
            if neighbour.cost < curr_solution.cost:
                curr_solution = neighbour
                break
        # if no direction of descent anymore
        if solution.cost >= curr_solution.cost:
            break
    return curr_solution

def _shake_simple(solution, neighbourhood_type, use_c):
    return Solution(random.choice(solution.get_neighbourhood(neighbourhood_type)), 
                    solution.problem,
                    use_c)

def shake(solution, neighbourhood_type, use_c, k=1):
    # For k > 1 shaking is significantly less eficient, because cost caluclations are needed.
    # That's why for k == 1 we still want to use "simple" shaking.
    if k == 1:
        return _shake_simple(solution, neighbourhood_type, use_c=use_c)

    neighbourhood = [Solution(neighbour, solution.problem, use_c)
                    for neighbour in solution.get_neighbourhood(neighbourhood_type)]
    random_k_sample = random.sample(neighbourhood, k)
    return get_best_solution(random_k_sample)

def basic_VNS(problem,
              initial_solution=None,
              diversification_param=0,
              initialization_method=get_initial_solution_robust,
              local_search=local_search_best_improvement,
              use_c=True,
              neighbourhood_types=NEIGHBOURHOOD_TYPES[:1],
              max_iter=MAX_ITER,
              precision=PRECISION,
              verbose=False):
    '''Basic VNS for UMApHMP.

    Parameters:
    problem (ProblemInstance)
    diversification_param (float): From [0, 1] interval. Describes the diversification of the search, \
        meaning 0 - most diversified, 1 - least diversified.
    precision (float): optimal solution precision
    verbose (bool)

    Notes:
    - Regarding diversification parameter:
    "In developing a more effective VNS, one must spend some time in checking how sensitive
    is the objective function to small change (shaking) of the solution.
    The trade-off between intensification and diversification of the search in VNS is balanced
    in a Shaking procedure. For some problem instances, a completely random jump in the kth neighborhood
    is too diversified. In such cases, some intensify shaking procedure is in order.
    For instance, a k-interchange neighbourhood may be reduced by repeating k times random add followed by best drop moves."
    Variable neighborhood search: Methods and applications, Pierre Hansen, Nenad Mladenovic and Jose A. Moreno Perez
    '''
    starttime = time.time()

    if diversification_param < 0 or diversification_param > 1:
        raise ValueError("diversification_param is expected to be from [0, 1] interval.")

    # get shaking intensity from diversification_param
    # the formula bellow maps [0, 1] into {1, 2... p(n-p)}; p(n-p) == len(neighbourghood)
    shake_param_k = int((problem.p*(problem.n - problem.p) - 1)*diversification_param + 1)

    # initialize solution
    if initial_solution == None:
        solution = Solution(initialization_method(problem.n, problem.p, problem.distances), problem, use_c=use_c)
    else:
        solution = initial_solution

    optimal_solution = solution
    if verbose == True:
        iter_range = tqdm(range(max_iter))
    else:
        iter_range = range(max_iter)
    for iteration in iter_range:
        i = 0
        while i < len(neighbourhood_types):
            # Shaking
            rand_solution = shake(solution, neighbourhood_types[i], use_c, shake_param_k)
            # Local search
            # let's do local search always just in the smalles neighbourhood, but shake in wider neighbourhoods
            local_min = local_search(rand_solution, neighbourhood_types[0], use_c)
            # Change neighbourhood
            if local_min.cost < solution.cost:
                solution = local_min
                # reset loop
                i = 0
            else:
                i += 1
        # Update optimal solution
        if solution.cost < optimal_solution.cost:
            optimal_solution = solution

        # if we know the optimal solution and we reached it we stop
        if problem.optimal_cost is not None \
            and abs(optimal_solution.cost - problem.optimal_cost) < precision:
            break

        # elapsed time as additional stopping criterion
        now = time.time()
        if now-starttime >= T_MAX:
            break

    return optimal_solution


def reduced_VNS(problem,
                initial_solution=None,
                diversification_param=0.3,
                initialization_method=get_initial_solution_robust,
                use_c=True,
                neighbourhood_types=NEIGHBOURHOOD_TYPES[:2],
                max_iter=MAX_ITER,
                precision=PRECISION,
                verbose=False):
    """ Reduced VNS.
    Here, unlike in other variants, max_iter is limit for max number of iterations between two improvenments.
    """
    starttime = time.time()

    if diversification_param < 0 or diversification_param > 1:
        raise ValueError("diversification_param is expected to be from [0, 1] interval.")

    # get shaking intensity from diversification_param
    # the formula bellow maps [0, 1] into {1, 2... p(n-p)}; p(n-p) == len(neighbourghood)
    shake_param_k = int((problem.p*(problem.n - problem.p) - 1)*diversification_param + 1)

     # initialize solution
    if initial_solution == None:
        solution = Solution(initialization_method(problem.n, problem.p, problem.distances), problem, use_c=use_c)
    else:
        solution = initial_solution

    optimal_solution = solution

    max_iters_without_improvement = 0
    iters_without_improvement = 0

    while max_iters_without_improvement < max_iter:
        i = 0
        while i < len(neighbourhood_types):
            # Shaking
            local_min = shake(solution, neighbourhood_types[i], use_c, shake_param_k)
            # Change neighbourhood
            if local_min.cost < solution.cost:
                solution = local_min
                # reset loop
                i = 0
            else:
                i += 1
        # Update optimal solution
        if solution.cost < optimal_solution.cost:
            optimal_solution = solution
            iters_without_improvement = 0
        else:
            iters_without_improvement += 1

        max_iters_without_improvement = max(iters_without_improvement, max_iters_without_improvement)
        # if we know the optimal solution and we reached it we stop
        if problem.optimal_cost is not None \
            and abs(optimal_solution.cost - problem.optimal_cost) < precision:
            break

        # elapsed time as additional stopping criterion
        now = time.time()
        if now-starttime >= T_MAX:
            break

    return optimal_solution


def deterministic_VNS(problem,
                        initial_solution=None,
                        initialization_method=get_initial_solution_robust,
                        local_search=local_search_best_improvement,
                        use_c=True,
                        neighbourhood_types=NEIGHBOURHOOD_TYPES[:1],
                        max_iter=MAX_ITER,
                        precision=PRECISION,
                        verbose=False):

    starttime = time.time()

    # initialize solution
    if initial_solution == None:
        solution = Solution(initialization_method(problem.n, problem.p, problem.distances), problem, use_c=use_c)
    else:
        solution = initial_solution

    optimal_solution = solution

    if verbose == True:
        iter_range = tqdm(range(max_iter))
    else:
        iter_range = range(max_iter)
    for iteration in iter_range:
        i = 0
        while i < len(neighbourhood_types):
            # local search
            local_min = local_search(solution, neighbourhood_types[i], use_c)
            # Change neighbourhood
            if local_min.cost < solution.cost:
                solution = local_min
                # reset loop
                i = 0
            else:
                i += 1
        # Update optimal solution
        if solution.cost < optimal_solution.cost:
            optimal_solution = solution

        # if we know the optimal solution and we reached it we stop
        if problem.optimal_cost is not None \
            and abs(optimal_solution.cost - problem.optimal_cost) < precision:
            break

        # elapsed time as additional stopping criterion
        now = time.time()
        if now-starttime >= T_MAX:
            break

    return optimal_solution

def general_VNS(problem,
              initial_solution=None,
              diversification_param=0,
              initialization_method=get_initial_solution_robust,
              local_search=local_search_best_improvement,
              use_c=True,
              neighbourhood_types=NEIGHBOURHOOD_TYPES[:1],
              max_iter=MAX_ITER,
              precision=PRECISION,
              verbose=False):
    starttime = time.time()

    if diversification_param < 0 or diversification_param > 1:
        raise ValueError("diversification_param is expected to be from [0, 1] interval.")

    # get shaking intensity from diversification_param
    # the formula bellow maps [0, 1] into {1, 2... p(n-p)}; p(n-p) == len(neighbourghood)
    shake_param_k = int((problem.p*(problem.n - problem.p) - 1)*diversification_param + 1)

    # initialize solution
    if initial_solution == None:
        solution = Solution(initialization_method(problem.n, problem.p, problem.distances), problem, use_c=use_c)
    else:
        solution = initial_solution

    optimal_solution = solution
    if verbose == True:
        iter_range = tqdm(range(max_iter))
    else:
        iter_range = range(max_iter)
    for iteration in iter_range:
        i = 0
        while i < len(neighbourhood_types):
            # Shaking
            rand_solution = shake(solution, neighbourhood_types[i], use_c, shake_param_k)
            # Local search
            # let's do local search always just in the smalles neighbourhood, but shake in wider neighbourhoods
            local_min = deterministic_VNS(problem, initial_solution=rand_solution)
            # Change neighbourhood
            if local_min.cost < solution.cost:
                solution = local_min
                # reset loop
                i = 0
            else:
                i += 1
        # Update optimal solution
        if solution.cost < optimal_solution.cost:
            optimal_solution = solution

        # if we know the optimal solution and we reached it we stop
        if problem.optimal_cost is not None \
            and abs(optimal_solution.cost - problem.optimal_cost) < precision:
            break

        # elapsed time as additional stopping criterion
        now = time.time()
        if now-starttime >= T_MAX:
            break

    return optimal_solution