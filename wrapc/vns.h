double normal_paths_calculation(int n, double** demand, double** cost_graph, int** predecessors);
int** floyd_warshall_c_impl(double** graph, int* hubs, int n, int p);
double get_solution_cost_c_impl(int* hubs, int n, int p, double** distances, double** demand, double alpha, double delta, double ksi);
