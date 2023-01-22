#!/usr/bin/env python
# coding: utf-8

# # Tutorial 2.2: Deriving the EOM of the pendulum covered during the lecture
# In this tutorial you will learn to derive the Equations of Motion of the pendulum as covered in the lecture. That is, a pendulum attached to the floor with a spring. We also consider an external horizontal force acting on the mass.

# ## Part 1: Kinematic equations
# Using the following rotation matrix R($\phi$): https://en.wikipedia.org/wiki/Rotation_matrix
# The zero angle points to the right for this matrix.
# 
# We first start by defining the variables. Assuming that x1, z1 and r1 do not change in time:

# In[1]:


import numpy as np
from sympy import *

var("t s x1 z1 r1")      
phi1 = Function("phi1")(t)

# Define the kinematic relations here
x2 = x1 + r1*cos(phi1)
z2 = z1 + r1*sin(phi1)


# The velocities can then be obtained using:

# In[2]:


# Compute/define the velocities here
xdot = diff(x2, t)
zdot = diff(z2, t)


# ## Part 2: Energy equations
# ### Kinetic energy:
# The only mass in this system is the point mass located at P2 with mass $M$

# In[3]:


var("M")
# Define the kinetic energy here (T)
T = 0.5*M*(xdot**2 + zdot**2)


# This expression can be simplified to:

# In[4]:


T = simplify(T)
T.evalf()


# ### Potential energy:
# Assuming again that $l0$, the original length of the spring, does not change in time, then the elongation of the spring $dl$ is given by:

# In[5]:


var("dx dz l0 x3 z3")
# Define the spring elongation here dl(l0,x2,z2,x3,z3)
dx = x3-x2
dz = z3-z2
dl = sqrt(dx**2 + dz**2) - l0
dl.evalf()


# The work done by the spring between P2 and P3 with stiffness $k$ and gravity.

# In[6]:


var("k g")
# Define the potential energy here (V)
V = 0.5*k*dl**2 + M*g*z2
V = simplify(V)
V.evalf()


# ### Work by external force
# The work done by an external force working in the horizontal direction:

# In[7]:


Fx = Function("Fx")(t)

# Define your external work here (W)
W = x2*Fx
W.evalf()


# ## Step 3: Construct the Lagrangian

# In[8]:


# Define your Lagrangian here (L)
L = T - V - W
L.evalf()


# ## Step 4: Obtaining the EoM
# 
# In order to obtain the EoMs we have to take derivatives w.r.t. $\phi_1$ and its velocity. 

# In[9]:


# Compute the EOM here
EOM_phi = diff( diff(L, diff(phi1, t)), t) - diff(L, phi1)
EOM_phi = simplify(EOM_phi)
EOM_phi.evalf()


# Now we isolate it for the acceleration

# In[10]:


# var("phi0 epsilon")
# psi = Function("psi")(t) # perturbation function

# tmp1 = symbols("tmp1")
# tmp2 = symbols("tmp2")
# EOM_psi2 = EOM_phi.evalf(subs={diff(phi1, (t, 2)): tmp2, phi1: tmp1})
# EOM_psi2 = EOM_psi2.evalf(subs={tmp2: diff(phi0 + epsilon*psi, (t, 2)), tmp1: phi0 + epsilon*psi})
# EOM_psi2.evalf()

# EOM_lin = series(EOM_psi2, epsilon, n=2)
# EOM_lin.evalf()

# EOM_lin = EOM_lin.removeO().evalf(subs={epsilon: 1})
# EOM_lin.evalf()

# # Isolate the acceleration
# EOM_lin_iso = solve(EOM_phi, diff(psi, (t, 2)))
# EOM_lin_iso[0].evalf()


# The EOM of an actual pendulum can be recovered by removing the spring and the external force. Note that this returns a cos because $\phi (t)=0$ was set to the right, rather than downward which is more conventional for a pendulum.

# In[11]:


# # Setting k=0 you should recover the pendulum equation found in tutorial 2_1
# a = EOM_lin_iso[0].evalf(subs={phi0: 3*pi/2, k: 0})


# ## Bonus: Obtaining the EOMs for nDOF system
# 
# If your system contains multiple DOFs, you will have to take the derivative of $L$ towards each of these separately, thereby obtaining the EOM for each DOF. Let's say you obtained the following two EOMs (the EOMs are just random things I entered):

# In[12]:


var("m c k J d q")
u = Function("u")(t)
w = Function("w")(t)

eom1 = m*cos(w)*diff(u, (t, 2)) + m*sin(u)*diff(w, (t, 2)) - c*diff(w, t) - c*diff(w, t)**2*u*ln(sqrt(t))**2*exp(u*w) - k*u
eom1.evalf()


# In[13]:


eom2 = J*w**2*diff(u, (t, 2)) + J*diff(w, t)*diff(w, (t, 2)) - d*diff(w, t) + q*u
eom2.evalf()


# You can collect both EOMs into matrix form using the following code:

# In[14]:


# linear_eq_to_matrix only accepts symbols, therefore we have to substitute the second derivative functions with symbols a1 and a2
a1, a2 = symbols("a1 a2") # acceler u and acceler w
eom1 = eom1.evalf(subs={diff(u, (t, 2)): a1, diff(w, (t, 2)): a2})
eom2 = eom2.evalf(subs={diff(u, (t, 2)): a1, diff(w, (t, 2)): a2})
MTRX = linear_eq_to_matrix([eom1, eom2], [a1, a2])


# Resulting in the following system of equations:

# In[15]:


M = MTRX[0]
M.evalf()


# And the following right-hand-side:

# In[16]:


F = MTRX[1]
F.evalf()

