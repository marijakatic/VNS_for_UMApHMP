import matplotlib.pyplot as plt 
import numpy as np
import string

def plot_solution(points, hubs, flow, title='', point_labels='numbers', plot_all_lines=False, verbose=0):
  fig, ax = plt.subplots(1, figsize=(6, 6))
  _plot_solution(ax, points,  hubs, flow, title, point_labels, plot_all_lines, verbose)
  plt.show()

def plot_two_solutions(points, hubs1, hubs2, flow1, flow2, title1='', title2='', point_labels='numbers', plot_all_lines=False, verbose=0):
  fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
  _plot_solution(ax1, points, hubs1, flow1, title1, point_labels, plot_all_lines, verbose)
  _plot_solution(ax2, points, hubs2, flow2, title2, point_labels, plot_all_lines, verbose)
  plt.show()

# Plot solution as a Graph in Cartesian coordinate system
def _plot_solution(ax, points, hubs, flow, title='', point_labels='numbers', plot_all_lines=False, verbose=0):
  n = len(points)
  # Get x and y coordinates
  xs = [point[0] for point in points]
  ys = [point[1] for point in points]

  # Assign colors to points (hub red, non-hub blue)
  color_map = {0: 'blue', 1: 'red'}
  point_colors = [color_map[hub] for hub in hubs]

  # Set point names
  if point_labels == 'numbers':
    point_labels = list(range(n))
  elif point_labels == 'alphabet':
    point_labels = list(string.ascii_uppercase)
  else: 
    raise ValueError("Unknown point_labels. Supported point_labels are: 'numbers', 'alphabet'.")

  # Plot points
  ax.scatter(xs, ys, c=point_colors)
  for i in range(n):
      ax.annotate(point_labels[i], points[i])

  # Draw lines connecting points
  for i in range(n):
    for k in range(n):
      # Points a and b
      point_a = points[i]
      point_b = points[k]
      # Get x and y coordinates
      endpoints_x = [point_a[0], point_b[0]]
      endpoints_y = [point_a[1], point_b[1]]

      if verbose == 1:
        # Centers of each line
        AB_x = (point_a[0] + point_b[0])/2
        AB_y = (point_a[1] + point_b[1])/2
        # Shift center a bit if we already had the line 
        if i > k:
          AB_y += 0.2
        # Write flows at the middle of the line
        ax.text(AB_x, AB_y, '{0}→{1}: {2}'.format(point_labels[i], point_labels[k], flow[i,k]))

      # line width to represent flow
      if verbose == 2:
        flow_sum = np.sum(flow)
        eta = 5 * flow_sum/np.max(flow)
        # eta = 50
        linewidth = max(abs(round((eta * flow[i,k]) / flow_sum, 4)), 0.0001)
      else:
        linewidth = 1.5

      if plot_all_lines:
        # Plot tiny line between every pair of points
        ax.plot(endpoints_x, endpoints_y, c='pink', ls='-', lw=1.5, alpha=0.5)

      # Ephasize connections
      if flow[i,k] > 0:
        ax.plot(endpoints_x, endpoints_y, c='grey', ls='--', lw=linewidth, alpha=0.5)
      # Ephasize connections between hubs
      if hubs[i] == 1 and hubs[k] == 1:
        ax.plot(endpoints_x, endpoints_y, c='black', ls='-', lw=linewidth)
      
  ax.set_title(title)
  # Draw major and minor grid lines
  ax.grid(which='both', color='yellow', linewidth=1, linestyle='-', alpha=0.2)
