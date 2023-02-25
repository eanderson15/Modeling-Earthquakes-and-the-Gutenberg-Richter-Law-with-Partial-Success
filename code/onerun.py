'''
Author: Eric Anderson
Program details from N. J. Giordano and H. Nakanishi, 2005, Computational Physics, Addison-Wesley, New York.
See references from above source for information on associated research papers.
Last edited: April 25, 2020
'''

from numpy import *
from vpython import *
from random import *
from time import *

# constants
m = 1.0
kp = 40.0
kc = 250.0
f_0 = 50.0
v_0 = 0.01

# changeable constants
# number of blocks
n = 25
# long timestep
dt1 = 0.03
# short timestep
dt2 = 0.003
# initial time
t_i = 0.0
# final time
t_f = 1000
# window size
window_w = 400
window_h = 400
# plus or minus factor to initial positions
rand = 0.001
# use uniformly distributed start or not
equil = False

# returns force of friction
def calculate_friction(f_opposed, v_i_t):
    # if block not moving
    if (v_i_t == 0.0):
        # if force on block less than equal to f_0
        if (f_opposed <= f_0):
            return -f_opposed
        # if force on block greater than f_0
        else:
            return sign(-f_opposed) * f_0
    else:
        # calculate kinetic friction
        return (-((f_0 * sign(v_i_t)) / (1.0 + abs(v_i_t))))

# return array of zeros
def initialize_positions_uniform(positions):
    return positions

# return array of positions adding or subtracting rand from each one
def initialize_positions_random(positions, rand):
    for i in range(n):
        switch_value = randint(0, 1)
        if (switch_value == 0):
            positions[i] = rand
        else:
            positions[i] = -rand
    return positions
 
# return value for force on block i       
def calculate_force(i):
    # if first block then no block to the left
    if (i == 0):
        f_opposed = kc * (positions[i + 1] - positions[i])
    # if last block then no block to the right
    elif (i == n - 1):
        f_opposed = kc * (positions[i - 1] - positions[i])
    else:
        f_opposed = kc * (positions[i + 1] + positions[i - 1] - 2.0 * positions[i])
    # add force from leaf spring
    f_opposed += kp * (v_0 * t_i - positions[i])
    # calculate and add friction
    ff = calculate_friction(f_opposed, velocities[i])
    return f_opposed + ff

# returns value for velocity of block i
def calculate_velocity(i, f_i):
    # calculate velocity of next step
    v_i = (f_i * dt) / m + velocities[i]
    # if velocity is not zero and it changes direction
    if ((velocities[i] != 0.0) and (sign(velocities[i]) != sign(v_i))):
        return 0.0
    return v_i

# returns value for position of block i
def calculate_position(i, v_i_t1):
    return (v_i_t1 * dt) + positions[i]

# plots the position and velocity for all blocks
def plots(): 
    # get the minimum and maximum position
    min_y_pos = float('inf')
    max_y_pos = float('-inf')
    for i in position_curves_positions:
        for j in i:
            if (j.y < min_y_pos):
                min_y_pos = j.y
            if (j.y > max_y_pos):
                max_y_pos = j.y
    # get the minimum and maximum velocit
    min_y_vel = float('inf')
    max_y_vel = float('-inf')
    for i in velocity_curves_velocities:
        for j in i:
            if (j.y < min_y_vel):
                min_y_vel = j.y
            if (j.y > max_y_vel):
                max_y_vel = j.y
    # compute factors for scaling
    xfactor_pos = window_w / t_f
    yfactor_pos = window_h / (max_y_pos - min_y_pos)
    xfactor_vel = window_w / t_f
    yfactor_vel = window_h / (max_y_vel - min_y_pos)
    # set up canvas and curves
    position_all = canvas(center = vector(t_i / 2.0, (max_y_pos - min_y_pos) / 2.0, 0.0), align = 'left', width = window_w, height = window_h, title = 'Position versus time', xtitle = 'Time', ytitle = 'Position')
    position_curves = []
    for i in range(n):
        position_curves.append(curve(color = vector(1.0, float(i) / n, 0.0)))
    # plot positions
    for i in range(n):
        for j in position_curves_positions[i]:
            position_curves[i].append(pos = vector(j.x * xfactor_pos, j.y * yfactor_pos, j.z))
    # set up canvas and curves
    velocity_all = canvas(center = vector(t_i / 2.0, (max_y_vel - min_y_vel) / 2.0, 0.0), align = 'left', width = window_w, height = window_h, title = 'Velocity versus time', xtitle = 'Time', ytitle = 'Velocity')
    velocity_curves = []
    for i in range(n):
        velocity_curves.append(curve(color = vector(1.0, float(i) / n, 0.0)))
    # plot velocities
    for i in range(n):
        for j in velocity_curves_velocities[i]:
            velocity_curves[i].append(pos = vector(j.x * xfactor_vel, j.y * yfactor_vel, j.z))

# get start time
start = time()

# initialized constants
# array of positions for each block
positions = resize(array([0.0]), n)
# array of velocities for each block
velocities = resize(array([0.0]), n)
# array of forces for each block
forces = resize(array([0.0]), n)

dt = dt1

# initialize positions
if (equil):
    positions = initialize_positions_uniform(positions)
else:
    positions = initialize_positions_random(positions, rand)

# graph and curve for position of block 10
position10 = graph(align = 'left', width = window_w, height = window_h, title = 'Position versus time', xtitle = 'Time', ytitle = 'Position')
g_position10 = gcurve(color=color.cyan)

# graph and curve for velocity of block 10
velocity10 = graph(align = 'left', width = window_w, height = window_h, title = 'Velocity versus time', xtitle = 'Time', ytitle = 'Position')
g_velocity10 = gcurve(color=color.cyan)

# arrays to store positions and velocities of all blocks
position_curves_positions = []
velocity_curves_velocities = []
for i in range(n):
    position_curves_positions.append([])
    velocity_curves_velocities.append([])

# keeps track of switching dt
first = False
# keeps track of dt having been switched until quake starts
second = False

# for this amount of time
while (t_i < t_f):
    # copy arrays so as not to overwrite them
    temp_positions = positions.copy()
    temp_velocities = velocities.copy()
    temp_forces = forces.copy()
    
    # keeps track of all velocities being zero
    zeros = True
    
    # iterate through blocks
    for i in range(n):
        f_i = calculate_force(i)
        v_i_t1 = calculate_velocity(i, f_i)
        
        # if block is moving
        if (v_i_t1 != 0.0):
            # if using long timestep
            if (dt == dt1):
                # go back
                t_i = t_i - dt
                # switch to short timestep
                dt = dt2
                first = True
                zeros = False
                break
            # marks quake start
            if (second):
                second = False
            zeros = False
            
        x_i_t1 = calculate_position(i, v_i_t1)
        
        # put calculated values in arrays
        temp_forces[i] = f_i
        temp_velocities[i] = v_i_t1
        temp_positions[i] = x_i_t1
    
    # if not switch (when go back in time)
    # update arrays and plot
    if (not first):
        g_position10.plot(pos = [t_i, positions[9]])
        g_velocity10.plot(pos = [t_i, velocities[9]])
        for i in range(n):
            position_curves_positions[i].append(vector(t_i, positions[i], i))
            velocity_curves_velocities[i].append(vector(t_i, velocities[i], i))
        forces = temp_forces
        velocities = temp_velocities
        positions = temp_positions
        
    # after switched
    else:
        first = False
        second = True
    
    # marks that quake is over
    if (zeros and not second):
        # switch to long timestep
        dt = dt1
        
    t_i += dt
    
plots()
print('Finish: ', time() - start)