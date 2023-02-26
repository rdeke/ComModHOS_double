#!/usr/bin/env python
# coding: utf-8

# # 1.1. Taylor series
# 
# Before going into more details on how to find a numerical approximation, let's start by refreshing some theory about **[Taylor series](https://en.wikipedia.org/wiki/Taylor_series)**. As one can find in the Wikipedia page:
# 
# > *the **Taylor series** of a function is an infinite sum of terms that are expressed in terms of the function's derivatives at a single point.*
# 
# That is $$ f(x)=\sum_{i=0}^\infty\frac{f^{(i)}(a)}{i!}(x-a)^i. $$
# 
# The series is exact for an arbitrary function as long as we include infinite terms in the sum. However, here we are interested on an approximation that includes only the first $r$ terms of the expansion: $$ f_r(x)=\sum_{i=0}^r\frac{f^{(i)}(a)}{i!}(x-a)^i\approx f(x). $$
# 
# Let's see how this work in a practical example. We consider here the function $f(x)=\sin(x)$ and we want to approximate this function knowing the function value and its derivatives at the point $a=0$. Using five terms in the expansion, $r=5$, we have that 
# $$\begin{align*}f_5(x) &= \sin(0) + x\cos(0) - \frac{x^2}{2}\sin(0) - \frac{x^3}{6}\cos(0) + \frac{x^4}{24}\sin(0) + \frac{x^5}{120}\cos(0)\\
#                       &= x - \frac{x^3}{6} + \frac{x^5}{120}.
# \end{align*} $$
# 
# We now can do a first coding exercise with this example and see how this approximation looks like for different values of $r$. To do so, we explicitly define the approximated functions.
# ```python
# import numpy as np
# 
# n = 200 # Gives smoother functions than standard n=50 linspace
# x = np.linspace(-2*np.pi,2*np.pi,n)
# a_0 = 0
# a = a_0*np.ones(len(x)) # Lower 
# 
# def f(x):
#   return np.sin(x)
# def f0(x):
#   return np.sin(a)
# def f1(x):
#   return f0(x) + (x-a)*np.cos(a)
# def f2(x):
#   return f1(x) - (x-a)**2/2 * np.sin(a)
# def f3(x):
#   return f2(x) - (x-a)**3/6 * np.cos(a)
# def f4(x):
#   return f3(x) + (x-a)**4/24 * np.sin(a)
# def f5(x):
#   return f4(x) + (x-a)**5/120 * np.cos(a)
# ```
# Let's plot all these functions:
# ```python
# import matplotlib.pyplot as plt
# %matplotlib inline
# 
# plt.plot(x,f(x),label="sin(x)")
# plt.plot(x,f0(x),label="f_0(x)",ls=":")
# plt.plot(x,f1(x),label="f_1(x)")
# plt.plot(x,f2(x),label="f_2(x)",ls=(0, (5, 10)))
# plt.plot(x,f3(x),label="f_3(x)")
# plt.plot(x,f4(x),label="f_4(x)",ls="dashed")
# plt.plot(x,f5(x),label="f_5(x)")
# 
# plt.xlabel("x")
# plt.ylabel("f(x)"),
# plt.xlim((-2*np.pi,2*np.pi))
# plt.ylim((-1.5,1.5))
# plt.legend();
# ```
# 
# \fig{../../images/Module1/TaylorSeries/1_2.png}
# 
# We see that as we increase the number of terms, the approximation gets closer to the analytical expression. Note also that for this particular example, there is no gain for $r=2$ and $r=4$ with respect to $r=1$ and $r=3$. This is caused by the fact that $\sin(a)=0$ at the approximation point $a=0$.
# 
# We can go further and evaluate and plot the error of the approximation for different values of $m$.
# $$ e_m(x)=|f(x)-f_m(x)|,\qquad r=1,...,5$$
# 
# ```python
# e0 = abs(f(x)-f0(x))
# e1 = abs(f(x)-f1(x))
# e2 = abs(f(x)-f2(x))
# e3 = abs(f(x)-f3(x))
# e4 = abs(f(x)-f4(x))
# e5 = abs(f(x)-f5(x))
# 
# plt.plot(x,np.log10(e0),label="$e_0(x)$",ls=":")
# plt.plot(x,np.log10(e1),label="$e_1(x)$")
# plt.plot(x,np.log10(e2),label="$e_2(x)$",ls=(0, (5, 10)))
# plt.plot(x,np.log10(e3),label="$e_3(x)$")
# plt.plot(x,np.log10(e4),label="$e_4(x)$",ls="dashed")
# plt.plot(x,np.log10(e5),label="$e_5(x)$")
# 
# plt.xlabel("x")
# plt.ylabel("log10(e(x))")
# plt.xlim((-2*np.pi,2*np.pi))
# plt.ylim((10**(-2), 2))
# plt.yscale("log")
# plt.legend();
# 
# plt.savefig("../../images/Module1/TaylorSeries/1_3.png")
# ```
# \fig{../../images/Module1/TaylorSeries/1_2.png}
# 
# As it is seen in the previous figure, close to the evaluation point $a$, the approximation error decreases as we increase the number of terms in the expansion.
# 
# :::{important} *Now we know...*
# * How to approximate a function using a **Taylor series**
# * That the approximation error is introduced when truncating the series to $r$ finite terms
# * How to compute and plot the Taylor series approximation and error of a simple function
# :::
# 
# ## Additional material
# In practise, one can use tools that facilitate the derivation of Taylor series of functions. Primarily these all use the same structure as shown above, however the different functions are created in a loop that more closely resembles the mathematical formulation using factorials. The symbolic basis of the used derivatives is most often done with the [`Sympy package`](https://docs.sympy.org/latest/modules/series/series.html).
# 
# A Similar tool exists for Julia, another scientific programming language that we will look into later in the unit. Then we use the package [`TaylorSeries.jl`](https://juliadiff.org/TaylorSeries.jl/stable/), which is a Julia package for Taylor expansions in one or more independent variables.
