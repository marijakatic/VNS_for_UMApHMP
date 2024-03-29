import numpy as np
import re

def parse_solutions(filepath):
    solutions = {}
    with open(filepath, mode='r', encoding='utf-8-sig') as file:
        # read lines ignoring empty lines and comments
        lines = [line for line in file.readlines() if line.strip() and not line.startswith('#')]
        for sol, obj, hubs in _chunks(lines, 3):
            match = re.match('Solution for n=(.*), p=(.*) :', sol)
            n = int(match.group(1))
            p = int(match.group(2))
            match = re.match('Objective : (.*)', obj)
            if match is not None:
                objective = float(match.group(1))
            else:
                objective = None
            match = re.match('Hubs : (.*)', hubs)
            if match is not None:
                hubs = _convert_list(match.group(1).split(', '), int)
            else:
                hubs = None
            solutions[n,p] = {"objective": objective, "hubs": hubs}
    return solutions

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

from docx import Document
import textract


def parse_CAB(filepath):
    with open(filepath, mode='r', encoding='utf-8-sig') as file:
        # read lines ignoring empty
        lines = [line for line in file.readlines() if line.strip()]

        # first line contains parameters of the problem
        parameters = lines[0].split()
        N, p = _convert_list(parameters[:2], int)
        alpha, delta, ksi = _convert_list(parameters[2:], float)
        
        # next N lines contain either C
        C = []
        for line in lines[1:N+1]:
            C.append(list(map(float, line.split())))

        # next N lines contain NxN matrix containing flows between points
        W = []
        for line in lines[N+1:]:
            W.append(list(map(float, line.split())))

        return (N, p, alpha, delta, ksi, np.array(C), np.array(W))


def parse_CAB_docx(filepath):
    text = textract.process(filepath).decode('utf-8-sig')
    # read lines ignoring empty
    lines = [line for line in text.split('\n\n') if line.strip()]

    # first line contains parameters of the problem
    parameters = lines[0].split()
    N, p = _convert_list(parameters[:2], int)
    alpha, delta, ksi = _convert_list(parameters[2:], float)

    # next N lines contain either C
    C = []
    for line in lines[1:N+1]:
        C.append(list(map(float, line.split())))

    # next N lines contain NxN matrix containing flows between points
    W = []
    for line in lines[N+1:]:
        W.append(list(map(float, line.split())))

    return (N, p, alpha, delta, ksi, np.array(C), np.array(W))

def get_comparison_table_file_name(experimentation_topic, number_of_problems, latest_commit_id):
    return  f'./output/basic_VNS_{experimentation_topic}_{number_of_problems}_{latest_commit_id}.csv'

def _convert_list(l, elem_type):
    return list(map(elem_type, l))

def _chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]