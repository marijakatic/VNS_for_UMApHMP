import numpy as np

def parse_input(filepath, dataset):
    if dataset == 'AP':
        return parse_AP(filepath)
    elif dataset == 'CAB':
        return parse_CAB(filepath)
    else:
        raise ValueError("Unknown dataset. Supported datasets: \
                         'AP' (\"Australia Post\" dataset), \
                         'CAB' (\"Civil Aeronautics Board\" datset)")
    
def parse_AP(filepath):
    with open(filepath, mode='r', encoding='utf-8-sig') as file:
        # read lines ignoring empty
        lines = [line for line in file.readlines() if line.strip()]

        # first line contains parameters of the problem
        N = int(lines[0])
        
        # next N lines contain N points
        points = []
        for line in lines[1:N+1]:
            points.append(tuple(map(float, line.split())))
        
        # next N lines contain NxN matrix containing flows between points
        W = []
        for line in lines[N+1:2*N+1]:
            W.append(list(map(float, line.split())))

        p = int(lines[2*N+1])
        ksi = float(lines[2*N+2])
        alpha = float(lines[2*N+3])
        delta = float(lines[2*N+4])

        return (N, p, alpha, delta, ksi, np.array(points), np.array(W))

def parse_CAB(filepath):
    with open(filepath, mode='r', encoding='utf-8-sig') as file:
        # read lines ignoring empty
        lines = [line for line in file.readlines() if line.strip()]

        # first line contains parameters of the problem
        parameters = lines[0].split()
        N, p = _convert_list(parameters[:2], int)
        alpha, delta, ksi = _convert_list(parameters[2:], float)
        
        # next N lines contain NxN matrix containing flows between points
        W = []
        for line in lines[1:N+1]:
            W.append(list(map(float, line.split())))

        # next N lines contain either W or C, idk :(
        C = []
        for line in lines[N+1:]:
            C.append(list(map(float, line.split())))

        return (N, p, alpha, delta, ksi, np.array(C), np.array(W))



def _convert_list(l, elem_type):
    return list(map(elem_type, l))
