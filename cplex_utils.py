from docplex.mp.model import Model
import datetime
import logging
import numpy as np
from utils import Solution

def get_flow_from_XYZ(n, X_allocated, Y_allocated, Z_allocated):
    """Calculates the flow matrix, the flow between pairs of nodes, given allocated variables from 
    the solution of Ernst and Krishnamoorthy formulation

    Parameters:
    n (int): number of nodes
    X_allocated (): 
    Y_allocated (): 
    Z_allocated (): 
    """
    flow = np.zeros((n,n))

    # Add up Zik 
    for i in range(n):
        for k in range(n):
            flow[i,k] += Z_allocated[i,k]

    # Add up Xilj
    for i in range(n):
        for l in range(n):
            for j in range(n):
                flow[l,j] += X_allocated[i,l,j]

    # Add up Yikl
    for i in range(n):
        for k in range(n):
            for l in range(n):
                flow[k,l] += Y_allocated[i,k,l]

    return flow

def solve_with_CPLEX(problem):
    M, X, Y, Z, H = get_umaphmp_model(problem.n, problem.p, problem.alpha, problem.delta, problem.ksi, problem.distances, problem.demand)
    solution = M.solve(log_output=True)
    return Solution(None, problem, False, solution.objective_value)

def get_umaphmp_model(n, p, alpha, delta, ksi, C, W, formulation='EK', verbose=True):
    """ Creates CPLEX model in a given formulation
    
    Parameters:
    N (int): number of nodes
    p (int): number of hubs
    alpha, delta, ksi (float): discounts
    C (numpy.ndarray): distance matrix
    W (numpy.ndarray): demand matrix

    Returns:
    docplex.mp.model.Model: Model to be solved with CPLEX solver
    """
    if formulation == 'EK':
        return ernst_krishnamoorthy(n, p, alpha, delta, ksi, C, W, verbose=verbose)
    elif formulation == 'C':
        return campbell(n, p, alpha, delta, ksi, C, W, verbose=verbose)
    else:
        raise ValueError("Unknown formulation. \
                          Possible values: 'EK' (Ernst and Krishnamoorthy formulation), 'C' (Campbell formulation).")
    

def ernst_krishnamoorthy(n, p, alpha, delta, ksi, C, W, model_name='UMApHMP, Ernst and Krishnamoorthy', verbose=True):
    """ Creates CPLEX model in Ernst and Krishnamoorthy formulation """
    if verbose:
        logging.basicConfig(format='[%(asctime)s] %(message)s', level = logging.INFO)
    else:
        logging.basicConfig(format='[%(asctime)s] %(message)s', level = logging.WARNING)

    
    # create model instance, with a name
    M = Model(model_name)
    logging.info("Created model M")

    # define X^{i}_{lj}
    X = {(i,l,j): M.continuous_var(name=f'X_{i}_{l}_{j}') for i in range(n) for l in range(n) for j in range(n)}
    logging.info("Defined variables X")
    # define Y^{i}_{kl}
    Y = {(i,k,l): M.continuous_var(name=f'Y_{i}_{k}_{l}') for i in range(n) for k in range(n) for l in range(n)} 
    logging.info("Defined variables Y")
    # define Zik
    Z ={(i,k): M.continuous_var(name=f'Z_{i}_{k}') for i in range(n) for k in range(n)}
    logging.info("Defined variables Z")
    # define Hk
    H = {k: M.binary_var(name=f'H_{k}') for k in range(n)} 
    logging.info("Defined variables H")

    # (2)
    M.add_constraint(M.sum(H[k] for k in range(n)) == p) 
    logging.info("Defined constraints (2)")
    # (3)
    for i in range(n):
        M.add_constraint(M.sum(Z[i,k] for k in range(n)) == M.sum(W[i,j] for j in range(n))) 
    logging.info("Defined constraints (3)")
    # (4)
    for i in range(n):
        for j in range(n):
            M.add_constraint(M.sum(X[i,l,j] for l in range(n)) == W[i,j]) 
    logging.info("Defined constraints (4)")
    # (5)
    for i in range(n):
        for k in range(n):
            M.add_constraint(M.sum(Y[i,k,l] for l in range(n)) + M.sum(X[i,k,j] for j in range(n)) - M.sum(Y[i,l,k] for l in range(n)) - Z[i,k] ==  0) 
    logging.info("Defined constraints (5)")
    # (6)
    for i in range(n):
        for k in range(n):
            M.add_constraint(Z[i,k] <= H[k]*M.sum(W[i,j] for j in range(n))) 
    logging.info("Defined constraints (6)")
    # (7)
    for l in range(n):
        for j in range(n):
            M.add_constraint(M.sum(X[i,l,j] for i in range(n)) <= H[l]*M.sum(W[i,j] for i in range(n))) 
    logging.info("Defined constraints (7)")
    # (8)
    # X, Y, Z >= 0 and Hk is binary are already met 

    # (1)
    M.minimize(M.sum((ksi*M.sum(C[i,k]*Z[i,k] for k in range(n)) +
                   alpha*M.sum(M.sum(C[k,l]*Y[i,k,l] for l in range(n)) for k in range(n)) +
                     delta*M.sum(M.sum(C[l,j]*X[i,l,j] for j in range(n)) for l in range(n))) for i in range(n)))
    logging.info("Defined constraints (1)")

    return (M, X, Y, Z, H)



def campbell(n, p, alpha, delta, ksi, C, W, model_name='UMApHMP, Campbell', verbose=True):
    if verbose:
        logging.basicConfig(format='[%(asctime)s] %(message)s', level = logging.INFO)
    else:
        logging.basicConfig(format='[%(asctime)s] %(message)s', level = logging.WARNING)

    # create model instance, with a name
    M = Model(model_name, log_output=verbose)
    logging.info("Created model M")

    # define Xijkm
    X = {(i,j,k,m): M.binary_var(name=f'X_{i}_{j}_{k}_{m}') for i in range(n) for j in range(n) for k in range(n) for m in range(n)}
    logging.info("Defined variables X")
    # define Hk (i.e. Xkk)
    H = {k: M.binary_var(name=f'H_{k}') for k in range(n)}
    logging.info("Defined variables H")

    # (4)
    M.add_constraint(M.sum(H[k] for k in range(n)) == p) 
    logging.info("Defined constraints (4)")
    # (5)
    # Xik is 0 or 1 - I am not sure why would we need it when it doesn't appear anymore
    # (10)
    # Xijkm >= 0 is default behavior
    # (14)
    for i in range(n):
        for j in range(n):
            M.add_constraint(M.sum(M.sum(X[i,j,k,m] for m in range(n)) for k in range(n)) == 1) 
    logging.info("Defined constraints (14)")
    # (15)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for m in range(n):
                    M.add_constraint(X[i,j,k,m] <= H[k])
    logging.info("Defined constraints (15)")
    # (16)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for m in range(n):
                    M.add_constraint(X[i,j,k,m] <= H[m]) 
    logging.info("Defined constraints (16)")

    # (7)
    M.minimize(M.sum(M.sum(M.sum(M.sum(W[i][j]*X[i,j,k,m]*(ksi*C[i,k] + delta*C[m,j] + alpha*C[k,m]) 
                                       for m in range(n)) for k in range(n)) for j in range(n)) for i in range(n))) 
    logging.info("Defined constraints (7)")

    return M, X, H