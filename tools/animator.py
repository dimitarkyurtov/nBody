import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import sys
import re

fileName = "nBody.txt"
# Read the file
if (len(sys. argv) >= 1):
    fileName = sys.argv[1]

print(sys.argv[1])     

with open(fileName, 'r') as f:
    lines = f.readlines()

# Get the range from the first line
x_range = list(map(float, lines[0].strip().split(',')))
print(x_range)

# Initialize the figure and axes
fig, ax = plt.subplots()
ax.set_xlim([-10000, 10000])
ax.set_ylim([-10000, 10000])

# Function to parse a line into time and points
def parse_line(line):
    time, points_str = line.strip().split(':')
    points = re.findall(r'\((.*?)\)', points_str)
    points = [tuple(map(float, point.split(','))) for point in points]
    return float(time), points

# Animation update function
def update(num, lines, scat):
    time, points = parse_line(lines[num])
    print(time, points)
    scat.set_offsets(points)
    return scat,

# Create the scatter plot
scat = ax.scatter([], [])

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=range(1, len(lines)), fargs=(lines, scat),
                              interval=1, blit=True)

plt.show()
