#!/usr/bin/env python
# coding: utf-8

# # Tutorial 1: ODE solvers
# In this tutorial you will learn to solve an MCK (mass-damper-spring) using a Python ODE solver. The MCK has 1 DOF and consequently the state veector contains 2 entries; displacement and velocity. 
# 
# $$ \bold{q} = \begin{bmatrix} u \\ \dot{u} \end{bmatrix}$$
# 
# The Equation of Motion (EoM) is given by: $$ m\ddot{u} = -ku -c\dot{u} $$

# ## Part 1: definition of inputs
# We start by defining the numerical values of all parameters:

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
k =         1   # [N/m] 
c =         0.1 # [Ns/m]
m =         1   # [kg = Ns^2/m]
#extra
a_f = 1
w_f = 3
def f(t):
    return a_f*np.sin(w_f*t)


# We want to solve the problem in the interval $ t \in [0,10] $, and get the solution with a resolution of $ \Delta t = 0.01$. Then, the vector of evaluation points will be defined as:

# In[2]:


# Time interval
t_0 =       0       # initial time [s]
t_f =       10      # final time [s]
dt =        0.01    # time step size [s]

steps =     int(t_f/dt)  # integer number of steps [-]
tspan =     np.linspace(t_0,t_f,steps)   # vector of evaluation points [s]


# The initial conditions for this example will be: $ u(0) = 1.0 $ and $ \dot{u}(0) = 0.0$.

# In[3]:


# Initial conditions
init_disp = 0.01 # m
init_velo = 0 # m/s


# ## Part 2: the ODE solver
# We will use solve_ivp (available via the SciPy package) to solve our MCK system. solve_ivp implements the 5th order Runge-Kutta scheme mentioned in the lectures. Checking help for solve_ivp shows how to call it:
# 
# `scipy.integrate.solve_ivp(fun, t_span, y0, t_eval, **options)`
# 
# OUTPUT:
# - T: Evaluation points
# - q: containts the state $q$ at each time in the vector T
# 
# INPUT:
# - `fun`: our ode function as explained during the lecture. It will compute the time derivatives of $q: \dot{q}=\mathcal{F} (q)$.
# - `t_span`: list of the first and last time step
# - `y0`: our initial state / conditions, $q(0)$.
# - `t_eval`: times at which the solution should be stored. In our case this is the variable `tspan`
# - `**options`: Options for the solvers. Here you can set things like error tolerances, maximum time step, event functions, etc. Check the SciPy docs for details. You will need to use this during the first assignment.  
# 
# Let's ignore the options for now and use the `solve_ivp` as our solver.

# ## Part 3: ODE function
# Now we need to create our `fun`. In the help documentation of `solve_ivp` we can check for a short description on what the function does. It is stated that the function f(t,q) determines the differential equations. The function is called as:
# 
# $ \dot{\bold{q}} = $`fun`$(t,\bold{q}) $
# 
# Here, the variable $\dot{\bold{q}}$ is the time derivative of our current state, $t$ is the current time and $\bold{q}$ is our current state. As the solver requires this interface, we have to create our `fun` accordingly or we will get answers that have no physical meaning!
# 
# Now we are faced with the problem that we need to have access to our parameters inside our `fun`. But our `fun` can only have the arguments $(t,\bold{q})$. To overcome this problem we will use an anonymous function (AF) to pass along the parameters to our `fun`. In Python an anonymous function is also known as a Lambda function.
# 
# `AF = lambda -,-,- : fun(-,-,-)`
# 
# where you have to replace the dashed marks. Now we can use this AF in our call to the solver:
# 
# `[T,q] = solve_ivp(AF,tspan,q_0)`
# 
# You can also directly declare the AF in the call to the solver if you prefer:
# 
# `[T,q] = solve_ivp(lambda -,-,- : fun(-,-,-),tspan,q_0)`
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Create a `fun` function that can receive the time, the state variable and the parameters as arguments. Implement the ODE function, $\mathcal{F}$, for the 1DOF MCK system such that $\dot{q}=\mathcal{F} (q)$.
# 
# *Hint*: Use the EoM and the fact that $\bold{q}(1) = u$ and $\bold{q}(2) = \dot{u}$.
# 
# -----------------------------------------------------------------------------------------------------
# 

# In[4]:


# Solve the problem of part 3 here
params = [k,c,m]
def fun(time, state): #time only relevant for time-dependent force
    
    k = params[0]
    c = params[1]
    m = params[2]
    
    DOF = int(len(state)/2)
    udotdot = np.zeros(DOF)
    for i in range(DOF):
        udotdot[i] = - ( k*state[2*i] + c*state[DOF+i] ) / m + f(time)
    
    state_dot = np.zeros(len(state))
    for i in range(DOF):
        state_dot[i] =  state[DOF+i]
        state_dot[DOF+i] = udotdot[i]
    
    return state_dot


# ## Part 4: initial state
# Next we need to create `q_0`. Note that the solver does not know / care what each entry in `q` represents. All the solver does is integrate things! You assign meaning to the entries in the state `q` when you define the initial conditions in `q_0`. This means that if you want the first entry in your state to be the displacement, `q_0[0]` should be set to `init_disp`. If you want the velocities to be the first entry, `q_0[0]` should be equal to the `init_velo`. It is up to you to decide this. 
# 
# !! IMPORTANT !!
# The `q` you receive in your `fun` will have the same meaning as your `q_0`. This means if you assigned `q_0[0]` to be the `init_disp`, `q_n[0]` will be the current displacement. So make sure your usage of `q_n` inside the `fun` is consistent with your definition of `q_0` as otherwise you will get bogus results
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Create your `q_0`
# 
# *Hint*: Straight forward! Just make sure that the indices are consistent with what you wrote in `fun`.
# 
# -----------------------------------------------------------------------------------------------------

# In[5]:


# Solve the problem of part 4 here
q_0 = [init_disp, init_velo]


# ## Part 5: Solve
# Once everything works the solver will return T and q. Each row in q corresponds to your state at that time-step. you can then plot your results with:

# In[6]:


# Solving system

sol = solve_ivp(fun=fun,t_span=[t_0,t_f],y0=q_0, t_eval=tspan)
print(sol)

# Plotting the solution
plt.plot(sol.t,sol.y[0],label="Displacement")
plt.plot(sol.t,sol.y[1],label="Velocity")
plt.xlabel('Time [s]')
plt.legend();


# In[7]:


def FE_solver1(fun,tspan,y0):
    
    q = np.zeros((2,len(tspan)))
    q[0][0] = y0[0]
    q[1][0] = y0[1]
    for i in range(1,len(tspan)):
        q_new = fun(tspan[i],[q[0][i-1],q[1][i-1]])
        dt = tspan[i] - tspan[i-1]
        q[0][i] = q[0][i-1] + q_new[0]*dt + 0.5*q_new[1]*dt**2
        q[1][i] = q_new[0] + q_new[1]*dt
        
    return q

# Plotting the solution
q = FE_solver1(fun,tspan,q_0)
plt.plot(tspan,q[0],label="Displacement")
plt.plot(tspan,q[1],label="Velocity")
plt.xlabel('Time [s]')
plt.legend();


# In[8]:


# Try ussing Forward-Euler

def FE_solver2(fun,tspan,y0):

    q = np.zeros((len(tspan),2))
    q[0] = y0
    for i,t in enumerate(tspan[1:]):
        q_dot_i = fun(t,q[i])
        q[i+1] = q[i] + np.multiply(q_dot_i,dt)
        
    return q

# Plotting the solution
q = FE_solver2(fun,tspan,q_0)
x = [q[i][0] for i in range(len(tspan))]
v = [q[i][1] for i in range(len(tspan))]
plt.plot(tspan,x,label="Displacement")
plt.plot(tspan,v,label="Velocity")
plt.xlabel('Time [s]')
plt.legend();


# In[9]:


w_n = np.sqrt(k/m)
zeta = c/(2*m*w_n)
w_d = w_n*np.sqrt(1-zeta**2)
A =  a_f * (m-k*w_f**2) / ( (m-k*w_f**2)**2 + (c*w_f)**2 )
B = -a_f * (c*w_f)      / ( (m-k*w_f**2)**2 + (c*w_f)**2 )
C_1 = (init_velo - A*w_f + zeta*w_n*(init_disp-B))/w_d
C_2 = init_disp - B

x_p = np.zeros(len(tspan))
x_h = np.zeros(len(tspan))
x = np.zeros(len(tspan))
v = np.zeros(len(tspan))
for i in range(len(tspan)):
    t = tspan[i]
    x_p[i] = A*np.sin(w_f*t) + B*(w_f*t)
    x_h[i] = np.exp(-zeta*w_n*t) * (C_1*np.sin(w_d*t) + C_2*np.cos(w_d*t))
    x[i] = x_p[i] + x_h[i]
for i in range(1,len(tspan)):
    v[i] = (x[i] - x[i-1])/(tspan[i] - tspan[i-1])
    
plt.plot(tspan,x,label="Analytical Displacement")
plt.plot(tspan,v,label="Velocity Displacement")
plt.xlabel('Time [s]')
plt.legend();


# In[10]:


# Error plotting (between 2 different numerical methods)
plt.plot(tspan,abs(x-sol.y[0]),label="Displacement")
plt.plot(tspan,abs(v-sol.y[1]),label="Velocity")
plt.xlabel('Time [s]')
plt.legend();

