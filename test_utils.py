import unittest
import utils

from os import listdir
from os.path import isfile

import ioutils


DATASET = 'AP'
INPUT_DIRECTORY = f"./data/{DATASET}/generated/"
SOLUTIONS_FILE = '/home/ubuntu/VNS_for_UMApHMP/data/AP/Solutions-UMApHMP.txt'

class UtilsTestCase(unittest.TestCase):

    def test_compare_get_solution_cost_and_get_solution_cost_fw(self):
        files = listdir(INPUT_DIRECTORY)
        files.sort()
        filter(isfile, files)
        print(files)

        for file in files:
            n, p, alpha, delta, ksi, points, demand = ioutils.parse_input(INPUT_DIRECTORY + file, DATASET)
            problem = utils.ProblemInstance(n, p, alpha, delta, ksi, points, demand)

            # solutions[n,p] = {"objective": objective, "hubs": hubs}
            solutions = ioutils.parse_solutions(SOLUTIONS_FILE)
            if (n, p) not in solutions or 'hubs' not in solutions[n, p]:
                continue
            hubs = solutions[n,p]['hubs']

            dijkstra_cost = round(utils.get_solution_cost(hubs, problem), 4)
            fw_cost = round(utils.get_solution_cost_fw(hubs, problem), 4)
            
            self.assertEqual(dijkstra_cost, fw_cost)


if __name__ == '__main__':
    unittest.main()