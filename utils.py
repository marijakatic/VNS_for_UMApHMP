from scipy.spatial import distance_matrix as _distance_matrix
import numpy as np
from graph_utils import dijkstra


def get_distance_matrix(nodes_coordinates):
    return _distance_matrix(nodes_coordinates, nodes_coordinates)

def get_discount_matrix(n, hubs, alpha, delta, ksi):
    discounts = np.ones((n,n))
    nodes = list(range(n))
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

def allocate_paths(n, hubs, distances, discounts):
    nodes = list(range(n))
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