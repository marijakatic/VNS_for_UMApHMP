from scipy.spatial import distance_matrix as _distance_matrix
import numpy as np
from graph_utils import dijkstra

from os.path import isfile

from ioutils import parse_input
from ioutils import parse_solutions
from plot_utils import plot_two_solutions

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
    

def allocate_paths(n, hubs, distances, discounts):
    nodes = get_nodes(n)
    cost_graph = distances*discounts
    for A in nodes:
        for B in nodes:
            # connection between non-hub nodes is not allowed
            if A != B and (A not in hubs) and (B not in hubs):
                cost_graph[A, B] = 0
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

def get_swap_neighbourhood(hubs, n):
    neighbourhood = []
    non_hubs = [node for node in range(n) if node not in hubs]
    for hub in hubs:
        for non_hub in non_hubs:
            # swap the hub with the non-hub
            neighbourhood.append([h for h in hubs if h != hub] + [non_hub])
    return neighbourhood

def bitmap(n, list):
    return [1 if i in list else 0 for i in range(n)]

def plot_comparison_with_optimal(input_directory, solutions_file, dataset, get_solution):
    '''
    Parameters:
    input_directory (str): directory with input files
    solution_file (str): file with solutions
    dataset (str): dataset label
    get_solution (method): method that returns a solution (as a list of hubs)
    '''
    solutions = parse_solutions(solutions_file)
    for n, p in solutions:
        # parse the instance n, p
        filepath = input_directory + f'{n}.{p}'
        if not isfile(filepath):
            continue
        n, p, alpha, delta, ksi, nodes_coordinates, demand = parse_input(filepath, dataset)
        distances = get_distance_matrix(nodes_coordinates)

        optimal_hubs = solutions[n,p]['hubs']
        optimal_discounts = get_discount_matrix(n, optimal_hubs, alpha, delta, ksi)
        optimal_paths = allocate_paths(n, optimal_hubs, distances, optimal_discounts)

        initial_hubs = get_solution(n, p, distances)
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