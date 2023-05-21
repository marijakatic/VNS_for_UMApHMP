from docplex.mp.model import Model
import datetime

def UMApHMP_f2(N, p, alpha, delta, ksi, C, W, model_name='UMApHMP_f2', verbose=True):
    # create model instance, with a name
    M = Model(model_name)
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Created model M")

    range_N = range(N)

    # define X^{i}_{lj}
    X = {(i,l,j): M.continuous_var(name='X_{0}_{1}_{2}'.format(i,l,j)) for i in range_N for l in range_N for j in range_N}
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined variable X")
    # define Y^{i}_{kl}
    Y = {(i,k,l): M.continuous_var(name='Y_{0}_{1}_{2}'.format(i,k,l)) for i in range_N for k in range_N for l in range_N} 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined variable Y")
    # define Zik
    Z ={(i,k): M.continuous_var(name='Z_{0}_{1}'.format(i,k)) for i in range_N for k in range_N}
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined variable Z")
    # define Hk
    H = {k: M.binary_var(name='H_{0}'.format(k)) for k in range_N} 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined variable H")

    # (2)
    M.add_constraint(M.sum(H[k] for k in range_N) == p) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (2)")
    # (3)
    for i in range_N:
        M.add_constraint(M.sum(Z[i,k] for k in range_N) == M.sum(W[i,j] for j in range_N)) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (3)")
    # (4)
    for i in range_N:
        for j in range_N:
            M.add_constraint(M.sum(X[i,l,j] for l in range_N) == W[i,j]) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (4)")
    # (5)
    for i in range_N:
        for k in range_N:
            M.add_constraint(M.sum(Y[i,k,l] for l in range_N) + M.sum(X[i,k,j] for j in range_N) - M.sum(Y[i,l,k] for l in range_N) - Z[i,k] ==  0) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (5)")
    # (6)
    for i in range_N:
        for k in range_N:
            M.add_constraint(Z[i,k] <= H[k]*M.sum(W[i,j] for j in range_N)) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (6)")
    # (7)
    for l in range_N:
        for j in range_N:
            M.add_constraint(M.sum(X[i,l,j] for i in range_N) <= H[l]*M.sum(W[i,j] for i in range_N)) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (7)")
    # (8)
    # X, Y, Z >= 0 and Hk is binary are already met 

    # (1)
    M.minimize(M.sum((ksi*M.sum(C[i,k]*Z[i,k] for k in range_N) +
                   alpha*M.sum(M.sum(C[k,l]*Y[i,k,l] for l in range_N) for k in range_N) +
                     delta*M.sum(M.sum(C[l,j]*X[i,l,j] for j in range_N) for l in range_N)) for i in range_N))
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (1)")

    return (M, X, Y, Z, H)



def UMApHMP_f1(N, p, alpha, delta, ksi, C, W, model_name='UMApHMP_f1', verbose=True):
    # create model instance, with a name
    M = Model(model_name)
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Created model M")

    range_N = range(N)

    # define Xijkm
    X = {(i,j,k,m): M.binary_var(name='X_{0}_{1}_{2}_{3}'.format(i,j,k,m)) for i in range_N for j in range_N for k in range_N for m in range_N}
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined variable X")
    # define Hk (i.e. Xkk)
    H = {k: M.binary_var(name='H_{0}'.format(k)) for k in range_N}
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined variable H")

    # (4)
    M.add_constraint(M.sum(H[k] for k in range_N) == p) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (4)")
    # (5)
    # Xik is 0 or 1 - I am not sure why would we need it when it doesn't appear anymore
    # (10)
    # Xijkm >= 0 is default behavior
    # (14)
    for i in range_N:
        for j in range_N:
            M.add_constraint(M.sum(M.sum(X[i,j,k,m] for m in range_N) for k in range_N) == 1) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (14)")
    # (15)
    for i in range_N:
        for j in range_N:
            for k in range_N:
                for m in range_N:
                    M.add_constraint(X[i,j,k,m] <= H[k])
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (15)")
    # (16)
    for i in range_N:
        for j in range_N:
            for k in range_N:
                for m in range_N:
                    M.add_constraint(X[i,j,k,m] <= H[m]) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (16)")

    # (7)
    M.minimize(M.sum(M.sum(M.sum(M.sum(W[i][j]*X[i,j,k,m]*(ksi*C[i,k] + delta*C[m,j] + alpha*C[k,m]) 
                                       for m in range_N) for k in range_N) for j in range_N) for i in range_N)) 
    if verbose:
        print("[" + str(datetime.datetime.now()) + "] Defined constraint (7)")

    return M, X, H
    






