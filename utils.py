from scipy.spatial import distance_matrix as _distance_matrix
import numpy as np
import pandas as pd

import sys
from os import listdir
from os.path import isfile
from os.path import basename
import time
from tqdm import tqdm
import git

from ioutils import parse_input
from ioutils import parse_solutions
from graph_utils import dijkstra
from graph_utils import floyd_warshall
from graph_utils import floyd_warshall_py
from graph_utils import NO_PATH_INDICATOR
from graph_utils import NO_EDGE_INDICATOR
from plot_utils import plot_two_solutions

from wrapc.wrapc import normal_paths_calulation_c
from wrapc.wrapc import floyd_warshall_c
from wrapc.wrapc import get_solution_cost_c

class ProblemInstance:
    def __init__(self, n, p, alpha, delta, ksi, node_coordinates, demand, distances=None, name=None, optimal_cost=None):
        self.n = n
        self.p = p
        self.alpha = alpha
        self.delta = delta
        self.ksi = ksi
        self.node_coordinates = node_coordinates
        self.demand = demand
        if distances is not None:
            self.distances = distances
        else:
            self.distances = get_distance_matrix(self.node_coordinates)
        self.optimal_cost = optimal_cost
        self.name = name

    def __str__(self):
        if self.name is None:
            return f"{self.n}.{self.p}"
        else:
            return self.name

    def __lt__(self, obj):
        return (self.n, self.p) < (obj.n, obj.p)

class Solution:
    def __init__(self, hubs, problem, use_c=True, cost=None):
        self.hubs = hubs
        self.problem = problem
        self.use_c = use_c
        if cost is not None:
            self.cost = cost
        else:
            self.cost = self._get_cost()
        # We'll caluclate the neighbourhood only when it's needed for the first time
        self.neighbourhood = None

    def __str__(self):
        return f"Solution(hubs={self.hubs})"

    def _get_cost(self):
        if self.use_c:
            return get_solution_cost_c(self.hubs, self.problem)
        else:
            return get_solution_cost_fw(self.hubs, self.problem, False)

    def get_neighbourhood(self, get_neighbourhood):
        if self.neighbourhood == None:
            self.neighbourhood = get_neighbourhood(self.hubs, self.problem.n)
        return self.neighbourhood


def get_nodes(n):
    return list(range(n))

def get_distance_matrix(nodes_coordinates):
    return _distance_matrix(nodes_coordinates, nodes_coordinates)

def get_discount_matrix(n, hubs, alpha, delta, ksi):
    discounts = np.ones((n,n))
    nodes = get_nodes(n)
    for A in nodes:
        for B in nodes:
            if (A in hubs) and (B in hubs):
                discounts[A, B] = alpha
            if (A not in hubs) and (B in hubs):
                discounts[A, B] = ksi
            if (A in hubs) and (B not in hubs):
                discounts[A, B] = delta
    return discounts

def get_total_cost(flow_matrix, distances, discounts):
    return sum(sum(flow_matrix*distances*discounts))

def get_flow_from_paths(n, paths, demand):
    flow = np.zeros((n,n))
    for path in paths:
        orig, dest = path[0], path[-1]
        for prev, cur in zip(path, path[1:]):
            flow[prev, cur] += demand[orig, dest]
    return flow

def get_flow_from_hubs(hubs, n, alpha, delta, ksi, distances, demand):
    discounts = get_discount_matrix(n, hubs, alpha, delta, ksi)
    paths = allocate_paths(n, hubs, distances, discounts)
    return get_flow_from_paths(n, paths, demand)
    
def _get_valid_cost_graph(n, hubs, distances, discounts):
    nodes = get_nodes(n)
    cost_graph = distances*discounts
    for A in nodes:
        for B in nodes:
            # connection between non-hub nodes is not allowed
            if A != B and (A not in hubs) and (B not in hubs):
                cost_graph[A, B] = NO_EDGE_INDICATOR
    return cost_graph

def allocate_paths(n, hubs, distances, discounts):
    nodes = get_nodes(n)
    # prepare graph
    cost_graph = _get_valid_cost_graph(n, hubs, distances, discounts)
    # calulate paths
    paths = []
    for A in nodes:
        _, paths_from_A = dijkstra(cost_graph, A)
        for B in nodes:
            if A != B:
                pathAB = paths_from_A[B]
                paths.append(pathAB)
    # by problem definition demand from a node to itself can be non-zero value
    # in that case, given demand should make a circuit through a hub and return to the node
    # although illogical at first glance, it could make sense in some setups
    # 
    # for example in post offices problem, if a non-hub node represents all residents of a city
    # a resident have to go to a post office (hub) to leave a mail, 
    # in order for it to be delivered to another resident of the same city
    more_paths = _orig_equals_dest_paths(nodes, hubs, cost_graph)
    return paths + more_paths

def get_solution_cost(hubs, problem):
    discounts = get_discount_matrix(problem.n, hubs, problem.alpha, problem.delta, problem.ksi)
    paths = allocate_paths(problem.n, hubs, problem.distances, discounts)
    flow = get_flow_from_paths(problem.n, paths, problem.demand)
    return get_total_cost(flow, problem.distances, discounts)

def _prepare(hubs, problem):
    nodes = get_nodes(problem.n)
    discounts = get_discount_matrix(problem.n, hubs, problem.alpha, problem.delta, problem.ksi)
    # prepare graph
    cost_graph = _get_valid_cost_graph(problem.n, hubs, problem.distances, discounts)
    # calculate cost
    return nodes, discounts, cost_graph

def _normal_paths_calulation(n, demand, cost_graph, predecessors):
    tmp_total_cost = 0
    for orig in range(n):
        for dest in range(n):
            left = predecessors[orig, dest]
            right = dest
            while True:
                if left == NO_PATH_INDICATOR:
                    break
                tmp_total_cost += demand[orig, dest] * cost_graph[left, right]
                if left == orig:
                    break
                right = left
                left = predecessors[orig, left]
    return tmp_total_cost

def _peculiar_paths_calulations(n, problem, hubs, cost_graph, discounts):
    tmp_total_cost = 0
    for node in range(n):
        if node in hubs:
            continue
        closest_hub = _closest(node, hubs, cost_graph)
        # adding cost of [node -> closest_hub -> node] path
        tmp_total_cost += problem.demand[node, node] * problem.distances[node, closest_hub] * discounts[node, closest_hub] + \
                      problem.demand[node, node] * problem.distances[closest_hub, node] * discounts[closest_hub, node]
    return tmp_total_cost

def get_solution_cost_fw(hubs, problem, use_c=True):
    total_cost = 0
    nodes, discounts, cost_graph = _prepare(hubs, problem)
    # calculate cost
    predecessors = floyd_warshall_c(cost_graph, hubs, problem.n, problem.p)

    if use_c:
        total_cost += normal_paths_calulation_c(problem.n, problem.demand, cost_graph, predecessors)
    else:
        total_cost += _normal_paths_calulation(problem.n, problem.demand, cost_graph, predecessors)

    # origin equals destination paths
    total_cost += _peculiar_paths_calulations(problem.n, problem, hubs, cost_graph, discounts)

    return total_cost

def get_swap_neighbourhood(hubs, n):
    neighbourhood = []
    non_hubs = [node for node in range(n) if node not in hubs]
    for hub in hubs:
        for non_hub in non_hubs:
            # swap the hub with the non-hub
            neighbourhood.append([h for h in hubs if h != hub] + [non_hub])
    return neighbourhood

def get_swap2_neighbourhood(hubs, n):
    neighbourhood = []
    non_hubs = [node for node in range(n) if node not in hubs]
    for hub1 in hubs:
        for hub2 in hubs:
            if hub1 == hub2:
                continue
            for non_hub1 in non_hubs:
                for non_hub2 in non_hubs:
                    if non_hub1 == non_hub2:
                        continue
                    # swap the hub with the non-hub
                    neighbourhood.append([h for h in hubs if h not in [hub1, hub2]] + [non_hub1, non_hub2])
    return neighbourhood

def get_swap3_neighbourhood(hubs, n):
    neighbourhood = []
    non_hubs = [node for node in range(n) if node not in hubs]
    for hub1 in hubs:
        for hub2 in hubs:
            for hub3 in hubs:
                if hub1 == hub2 or hub1 == hub3 or hub2 == hub3:
                    continue
                for non_hub1 in non_hubs:
                    for non_hub2 in non_hubs:
                        for non_hub3 in non_hubs:
                            if non_hub1 == non_hub2 or non_hub1 == non_hub2 or non_hub2 == non_hub3:
                                continue
                            # swap the hub with the non-hub
                            neighbourhood.append([h for h in hubs if h not in [hub1, hub2, hub3]] + [non_hub1, non_hub2, non_hub3])
    return neighbourhood

def bitmap(n, list):
    return [1 if i in list else 0 for i in range(n)]

def plot_comparison_with_optimal(input_directory, solutions_file, dataset, get_solution, num_of_problems):
    '''
    Parameters:
    input_directory (str): directory with input files
    solution_file (str): file with solutions
    dataset (str): dataset label
    get_solution (method): method that returns a solution (as a list of hubs)
    '''
    solutions = parse_solutions(solutions_file)
    problems = list(solutions.keys())[:num_of_problems]
    for n, p in problems:
        # parse the instance n, p
        filepath = input_directory + f'{n}.{p}'
        if not isfile(filepath):
            continue
        n, p, alpha, delta, ksi, nodes_coordinates, demand = parse_input(filepath, dataset)
        distances = get_distance_matrix(nodes_coordinates)

        optimal_hubs = solutions[n,p]['hubs']
        if optimal_hubs is None:
            continue
        optimal_discounts = get_discount_matrix(n, optimal_hubs, alpha, delta, ksi)
        optimal_paths = allocate_paths(n, optimal_hubs, distances, optimal_discounts)

        initial_hubs = get_solution(n, p, distances, nodes_coordinates)
        initial_discounts = get_discount_matrix(n, initial_hubs, alpha, delta, ksi)
        initial_paths = allocate_paths(n, initial_hubs, distances, initial_discounts)
        plot_two_solutions(nodes_coordinates,
                        bitmap(n, optimal_hubs),
                        bitmap(n, initial_hubs),
                        get_flow_from_paths(n, optimal_paths, demand),
                        get_flow_from_paths(n, initial_paths, demand),
                        title1=f"Optimal solution for n={n}, p={p}",
                        title2=f"Initial solution for n={n}, p={p}",
                        verbose=2)

def get_comparison_table(list_of_methods,
                         method_names,
                         dataset,
                         test_data_directory,
                         solutions_file=None,
                         from_problem=0,
                         to_problems=45,
                         initializations=[None],
                         display_solution=False,
                         display_deviation=True,
                         display_time=True,
                         display_nan=True):
    problems = []
    solutions = {}
    if solutions_file != None:
        solutions = parse_solutions(solutions_file)
    files = listdir(test_data_directory)
    for file in files:
        nodes = None
        distances = None
        if dataset == 'AP':
            n, p, alpha, delta, ksi, nodes, demand = parse_input(test_data_directory + file, dataset)
        elif dataset == 'CAB':
            n, p, alpha, delta, ksi, distances, demand = parse_input(test_data_directory + file, dataset)
        if (n,p) in solutions:
            optimal_cost = solutions[(n,p)]['objective']
        else:
            # todo fix this hack
            optimal_cost = None
        if dataset == 'AP':
            problem_name = None
        elif dataset == 'CAB':
            problem_name = basename(file)
        problems.append(ProblemInstance(n, p, alpha, delta, ksi, node_coordinates=nodes, distances=distances, demand=demand, optimal_cost=optimal_cost, name=problem_name))
    problems = sorted(problems)[from_problem:to_problems]

    if display_solution:
        columns = ['optimal solution']
    else:
        columns = []
    for i in range(len(list_of_methods)):
        method = list_of_methods[i]
        method_name = method_names[i]
        for initialization in initializations:
            if initialization != None and method_name != "cplex":
                if display_solution:
                    columns.append(method_name + " " + initialization.__name__ + " - solution")
                if display_deviation:
                    columns.append(method_name + " " + initialization.__name__ + " - deviation (%)")
                if display_time:
                    columns.append(method_name + " " + initialization.__name__ + " - time (s)")
            else:
                if display_solution:
                    columns.append(method_name + " - solution")
                if display_deviation:
                    columns.append(method_name + " - deviation (%)")
                if display_time:
                    columns.append(method_name + " - time (s)")

    data = {}
    for problem in tqdm(problems):
        if not display_nan and problem.optimal_cost == None:
            continue
        if display_solution:
            data[str(problem)] = [problem.optimal_cost]
        else:
            data[str(problem)] = []
        for method in list_of_methods:
            for initialization in initializations:
                if initialization != None and method.__name__ != "cplex":
                    initial_solution = Solution(initialization(problem.n, problem.p, problem.distances, problem.node_coordinates), problem)
                    starttime = time.time()
                    solution = method(problem, initial_solution=initial_solution)
                    endtime = time.time()
                else:
                    starttime = time.time()
                    solution = method(problem)
                    endtime = time.time()
                solution.cost = solution.cost/1000
                if display_solution:
                    data[str(problem)].append(round(solution.cost, 2))
                if problem.optimal_cost == None:
                    deviation = None
                else:
                    deviation = 100 * abs(problem.optimal_cost - solution.cost) / problem.optimal_cost
                    deviation = round(deviation, 2)

                if display_deviation:
                    data[str(problem)].append(deviation)
                if display_time:
                    data[str(problem)].append(round(endtime-starttime, 2))
    df = pd.DataFrame.from_dict(data, orient='index', columns=columns)
    df.loc['mean'] = df.mean()
    return df

def get_latest_commit_id():
    repo = git.Repo(search_parent_directories=True)
    return repo.head.object.hexsha

def _orig_equals_dest_paths(nodes, hubs, cost_graph):
    peculiar_paths = []
    for node in nodes:
        if node in hubs:
            continue
        closest_hub = _closest(node, hubs, cost_graph)
        peculiar_paths.append([node, closest_hub, node])
    return peculiar_paths

def _closest(node, hubs, cost_graph):
    hub_dist = {hub: cost_graph[node, hub] for hub in hubs}
    return min(hub_dist, key=hub_dist.get)



NEIGHBOURHOOD_TYPES = [get_swap_neighbourhood, get_swap2_neighbourhood, get_swap3_neighbourhood]
