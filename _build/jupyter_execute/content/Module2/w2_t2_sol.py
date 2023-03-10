#!/usr/bin/env python
# coding: utf-8

# # Solution 2.2: Deriving the EOM of a pendulum
# In this tutorial you will learn to derive the Equations of Motion of the pendulum as covered in the lecture. That is, a pendulum attached to the floor with a spring. We also consider an external horizontal force acting on the mass.
# 
# <center><img src="../../images/Module2/w2_t2_representation.png" width="200" /></center>

# ## Part 1: Kinematic equations
# Using the following rotation matrix R($\phi$): https://en.wikipedia.org/wiki/Rotation_matrix the zero angle points to the right for this matrix.
# 
# We first start by defining the variables. Assuming that x1, z1 and r1 do not change in time:

# In[1]:


import numpy as np
from sympy import *

var("t s x1 z1 r1")      
phi1 = Function("phi1")(t)


# -----------------------------------------------------------------------------------------------------
# **Problem**: Define the kinematic relations
# 
# *Hint*: How do you write the position of the pendulum (x2, z2) as a function of x1, z1, r1 and phi1?
# 
# ---------------------------------------------------------------------------------------------------

# In[2]:


# Define the kinematic relations here
x2 = x1 + r1*cos(phi1)
z2 = z1 + r1*sin(phi1)


# The velocities can then be obtained using:
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Define the velocities
# 
# *Hint*: Use derivatives.
# 
# ---------------------------------------------------------------------------------------------------

# In[3]:


# Compute/define the velocities here
xdot = diff(x2, t)
zdot = diff(z2, t)


# ## Part 2: Energy equations
# ### Kinetic energy:
# The only mass in this system is the point mass located at P2 with mass $M$
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Define the kinetic energy
# 
# ---------------------------------------------------------------------------------------------------

# In[4]:


var("M")
# Define the kinetic energy here (T)
T = 0.5*M*(xdot**2 + zdot**2)


# This expression can be simplified to:

# In[5]:


T = simplify(T)
T.evalf()


# ### Potential energy:
# Assuming again that $l0$, the original length of the spring, does not change in time, then the elongation of the spring $dl$ is given by:
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Define the spring elongation dl
# 
# ---------------------------------------------------------------------------------------------------

# In[6]:


var("dx dz l0 x3 z3")
# Define the spring elongation here dl(l0,x2,z2,x3,z3)
dx = x3-x2
dz = z3-z2
dl = sqrt(dx**2 + dz**2) - l0
dl.evalf()


# The work done by the spring between P2 and P3 with stiffness $k$ and gravity.
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Define the potential energy V.
# 
# ---------------------------------------------------------------------------------------------------

# In[7]:


var("k g")
# Define the potential energy here (V)
V = 0.5*k*dl**2 + M*g*z2
V = simplify(V)
V.evalf()


# ### Work by external force
# The work done by an external force working in the horizontal direction:
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Define the external work.
# 
# ---------------------------------------------------------------------------------------------------

# In[8]:


Fx = Function("Fx")(t)

# Define your external work here (W)
W = x2*Fx
W.evalf()


# ## Step 3: Construct the Lagrangian
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: Use T, V and W to find the lagrangian L. Simplify it using L.evalf()
# 
# ---------------------------------------------------------------------------------------------------

# In[9]:


# Define your Lagrangian here (L)
L = T - V - W
L.evalf()


# ## Step 4: Obtaining the EoM
# 
# In order to obtain the EoMs we have to take derivatives w.r.t. $\phi_1$ and its velocity. 
# 
# -----------------------------------------------------------------------------------------------------
# **Problem**: USe the Lagrangian expression to obtain the equation of motion.
# 
# *Hint*: Put all terms on the left hand side
# 
# ---------------------------------------------------------------------------------------------------

# In[10]:


# Compute the EOM here
EOM_phi = diff( diff(L, diff(phi1, t)), t) - diff(L, phi1)
EOM_phi = simplify(EOM_phi)
EOM_phi.evalf()


# Now we isolate it for the acceleration.
# 
# These cells can in practice take quite some time to execute. For this illustration the added spring stiffness k is set to 0, but it can take any value.

# In[11]:


var("phi0 epsilon")
psi = Function("psi")(t) # perturbation function

tmp1 = symbols("tmp1") # the acceleration, to be linearized as phi0+epsilon*psi
tmp2 = symbols("tmp2") # the displacement
EOM_psi2 = EOM_phi.evalf(subs={diff(phi1, (t, 2)): tmp2, phi1: tmp1})
EOM_psi2 = EOM_psi2.evalf(subs={tmp2: diff(phi0 + epsilon*psi, (t, 2)), tmp1: phi0 + epsilon*psi})
EOM_psi2.evalf()


# In[12]:


# Linearize by doing series expansion up to order 2 (see last term of output)
EOM_lin = series(EOM_psi2, epsilon, n=2)
EOM_lin.evalf()


# In[13]:


EOM_lin2 = EOM_lin.evalf(subs={k: 0}) # For faster execution
EOM_lin2.evalf()


# In[14]:


# Simplify and say that the "epsilon" was just a temorary value for perturbations
EOM_lin3 = EOM_lin2.removeO().evalf(subs={epsilon: 1})
EOM_lin3.evalf()


# In[15]:


# Isolate the acceleration
EOM_lin_iso = solve(EOM_lin3, diff(psi, (t, 2)))
EOM_lin_iso[0].evalf()


# The EOM of an actual pendulum can be recovered by removing the spring and the external force. Note that this returns a cos because $\phi (t)=0$ was set to the right, rather than downward which is more conventional for a pendulum.

# In[16]:


# Setting k=0 (done before) and phi0 straight down, you should recover the pendulum equation found in tutorial 2_1
EOM_lin_iso2 = EOM_lin_iso[0].evalf(subs={phi0: 3*pi/2, k: 0})
print(EOM_lin_iso2)


# In the code above the high length of the expressions is a consequence of machine precision. This makes the execution times a multitude larger as well. Depending on the application it becomes more and more usefulf to simplify as early as possible in hte process. addinh $k: 0$ to one of the earlier dictionaries for example does the trick here, but it can only be applied if the goal was to not include an additional spring in the first place.

# ## Bonus: Obtaining the EOMs for nDOF system
# 
# If your system contains multiple DOFs, you will have to take the derivative of $L$ towards each of these separately, thereby obtaining the EOM for each DOF. Let's say you obtained the following two EOMs (the EOMs are just random things I entered):

# In[17]:


var("m c k J d q")
u = Function("u")(t)
w = Function("w")(t)

eom1 = m*cos(w)*diff(u, (t, 2)) + m*sin(u)*diff(w, (t, 2)) - c*diff(w, t) - c*diff(w, t)**2*u*ln(sqrt(t))**2*exp(u*w) - k*u
eom1.evalf()


# In[18]:


eom2 = J*w**2*diff(u, (t, 2)) + J*diff(w, t)*diff(w, (t, 2)) - d*diff(w, t) + q*u
eom2.evalf()


# You can collect both EOMs into matrix form using the following code:

# In[19]:


# linear_eq_to_matrix only accepts symbols, therefore we have to substitute the second derivative functions with symbols a1 and a2
a1, a2 = symbols("a1 a2") # acceler u and acceler w
eom1 = eom1.evalf(subs={diff(u, (t, 2)): a1, diff(w, (t, 2)): a2})
eom2 = eom2.evalf(subs={diff(u, (t, 2)): a1, diff(w, (t, 2)): a2})
MTRX = linear_eq_to_matrix([eom1, eom2], [a1, a2])


# Resulting in the following system of equations:

# In[20]:


M = MTRX[0]
M.evalf()


# And the following right-hand-side:

# In[21]:


F = MTRX[1]
F.evalf()

