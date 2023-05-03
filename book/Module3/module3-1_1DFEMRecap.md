---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# 4.1. 1D FEM - a review

From earlier courses the general outline of 1D FEM should be clear. However, here a quick recap is goven, also to align any differences in notations or definitions we may use. 

The general goal of the FEM described here is to convert PDE's (Partial Differential Equations) into a form which we, or rather a computer, can solve. The reason for doing this is because it allows a PDE with no analytical solution for the full domain, to be approximated by solving it on smaller 'elements'. Because of this it becomes possible to solve just about any PDE using irregular bounds, different materials, dynamic effects, jumps or complex changes in geometry, ...

The general scheme for FEM problems is as follows:
1. Discretize the domain and get all elements, nodes, and their properties
2. Make piece-wise polynomials as approximations for the real solution
3. Define the weak form
4. Assemble the matrices to obtain the equations fpr the full solution
5. Sum interpolations for the solution

We will explain this for the case of a prismatic rod that is loaded in its axial direction. Here prismatic means that there no geometry or property changes along its length. The EOM that describes this phenomenon is 

$$ \rho A \frac{\partial^2 u(x,t)}{\partial t^2} - EA \frac{\partial^2u(x,t)}{\partial x^2} = q(x) $$

Additionally:
- We take the rod fixed at $x=0$, that means that$u(0,t) = 0 \forall t $. 
- The rod is at rest initially, $u(x,0) = 0$ and $\dot{u}(x,0) = 0$.
- A dynamic point load $P(t)$ is applied at the far end of the beam at $x=L$.

## 4.1.1 Step 1: Discretize the domain
For 1D examples discretization is very straightforward. The domain gets cut into $N$ pieces calles elements. Although not obligatory, most often each element has an equal length.

For a 1D case $N$ elements results in $N+1$ nodes, the boundaries of all elements.

## 4.1.2 Step 2: Define the shape functions
Given that the original FEM is second order, the shape functions must be at least order 1, so in this case linear polynomial functions. Other examples, as shown later, may require higher order polynomials.

We define the shape function $N_i^k$ as:
$$N_i^k(x) = a + bx$$
Here $i$ indicates that this is the shape function for the i'th element. The letter $k$ indicates that this is a function at 'element level'.

Shape functions are defined such that they have a value of 1 at node $i$. So mathematically this comes to
$$ N_i^k(x) = \begin{cases} 
          p(x) & x \in \Omega_k \forall k \in S_i \\
          1 & x = x_i \\
          0 & otherwise x 
       \end{cases}
$$
Here $p(x)$ just indicates a function defining values between 0 and 1. $\Omega$ itself denots the domain, with $\Omega_k$ the portion of the domain that lies around node $i$, so from node $i-1$ to node $i+1$.
In order to find the actual value of the displacement $u$ at any location we can use these shape functions to interpolate the the displacement from the known displacement in the two neigbouring nodes, so through:
$$ u(x) = \sum_{i=1}^{n_n} N_i^k(x) u_i$$
with $n_n$ being the number of nodes, so $N+1$.

 When looking at a single element bound by nodes $i$ and $i+1$ we thus get 2 equations from which we can solve the 2 unknowns:

$$ N_i^k  = a + bx 
            \Rightarrow
            \begin{cases} 
                N_i(x_i) = 1 \\
                N_i(x_{i+1}) = 0  \\ 
            \end{cases}
            \Rightarrow
            \begin{cases} 
                a + b x_i = 1 \\
                a + b x_{i+1} = 0  \\ 
            \end{cases}
            \Rightarrow
            \begin{cases} 
                a  = \frac{x_{i+1}}{x_{i+1}-x_{i}} \\
                b = \frac{1}{x_{i+1}-x_{i}} \\
              \end{cases}
$$

For element $i$ we need nodes $i$ and $i+1$. Mostly we write the shape functions with $h$ being the distance between the two nodes $x_{i+1}-x_i$$ so this gives

$$
\begin{cases} 
  N_i^k(x) & = \frac{x_{i+1} - x}{h} \\
  N_{i+1}^k(x) & = \frac{x-x_i}{h} \\
\end{cases}
$$

For shape functions of higher orders the way of obtaining $a, b, ...$ remains the same. However it can be practical to convert it into matrix form. This makes it easier when the number of unknowns grows

$$
\begin{bmatrix}
  1 & x_i \\
  1 & x_{i+1}
\end{bmatrix}
\begin{bmatrix}
  a \\
  b
\end{bmatrix}
=
\begin{bmatrix}
  1 \\
  0
\end{bmatrix}
$$

## 4.3 Define the weak form
The FEM algorithm relies on finding a balance through what is effectively the principle of virtual work. All nodes are displaced by $v$, which is a random *test* value. That gives the weak form. $v$ can in theory take any value we want, and then the EOM should be balanced.

We take the left hand side (LHS) of the equation, which is essentially a function of the displacement $u$ of each node as $L(u)$. The right hand side (RHS) we write generally as $f$:
$$ L(u) = f $$
Important to note here is that the RHS contains no dependencies on the nodal displacements $u$
We impose a virtual displacement $v$, and integrate both sides ofer the domain. Given we're at element level still, this equates to 
$$ \int_{\Omega_k} L(u) \cdot v d\Omega_k = \int_{\Omega_k} f \cdot v d\Omega_k $$
The LHS of this is often shortened to $a_k(u,v)$ and the RHS to $b_k(v)$.

To solve this equation we need 
$$ a_k(u,v)-b_k(v) = 0 \hspace{10pt} \forall v $$
The FEM relies then on the fact that $ u(x) = \sum_{i=1}^{n_n} N_i^k(x) u_i$

$$ a(u,v) - b(v) = 
\sum_{k=1}^{n_e} \left(a_k(u,v)-b_k(v)\right) = 
0 $$
Not that $n_e$ is the number of elements $N$.

Given then that $v$ can be any value, we take $v = N_i$ for $i$ in 1, ..., $n_n$. 

In case of the EOM described above this would look as follows:

$$

\int_{\Omega_k} m \ddot{u}(x) v(x) d\Omega_k - 
\int_{\Omega_k} EA u''(x) v(x) d\Omega_k =

\int_{\Omega_k} q(x) v(x) d\Omega_k

$$

We convert the second order derivative $u''$ through partial integration into a first order one. This makes the calculation simpler. This is the reason why we needed first order shape functions, since we reduce the derivative order of the original EOM to one. Partial integraton has the general form of 
$ \int_{x_i}^{x_{i+1}} f'' g d\Omega = -\int_{x_i}^{x_{i+1}} f' g' d\Omega+ f' g \vert_{x_i}^{x_{i+1}}
$. So this gives

$$

\int_{\Omega_k} m \ddot{u}(x) v(x) d\Omega_k - 
\int_{\Omega_k} EA u'(x) v'(x) d\Omega_k - 
EA u'(x) v(x) \vert_{x_i}^{x_{i+1}}
= \int_{\Omega_k} q(x) v(x) d\Omega_k

$$

Some interesting points to note here:
- The first integral contains the mass contribution
- The second integral contains the inertia contributions
- The third term describes the internal forces. 
- The RHS shos the external forces. Keen readers may already have realised that this value will be 0 for all nodes except the end-node (in our example)

As explained earlier we fill in $N_i$ for $v$, and then apply $u(x) = \sum_{i=1}^{n_n} N_i^k(x) u_i$, which brings us to

$$

\sum_{j=i,i+1}
\left(
\int_{\Omega_k} m N_j(x) N_i(x) d\Omega_k
\right) \ddot u_j + 
\sum_{j=i,i+1}
\left(\int_{\Omega_k} EA N_j'(x) N_i'(x) d\Omega_k - 
EA N_j'(x) N_i(x) \vert_{x_i}^{x_{i+1}}
\right) u_j
= \int_{\Omega_k} q(x) v(x) d\Omega_k

$$

This is a single equation with 2 unknowns, namely $u_i$ and $u_{i+1}$. However, this is only the element level. The neighbouring element $i+1$ relies on nodes $i+1$ and $i+2$. This makes a total system of $N$ equations (the number of elements) with $N+1$ unknowns (the number of nodes). Given that we know the displacement of one of the nodes already, through the boundary conditions, we can now assemble and then solve the system. 

The equation above can be rewritten into matrix form. With their dependency on $i$ and $j$ we write
$$ M_{ij}^k \ddot u_j + K_{ij}^k u_j = Q_i^k + S_i^k$$

Here $S_i^k = EA N_j'(x) N_i(x) \vert_{x_i}^{x_{i+1}}$. This term is brought to the RHS given that it does not actually depend on u. The term denotes the internal forces of the full system, or on alemental level the external forces on the element.

When we finally assemble the elemental system the full matrix equation then becomes
$$ \mathbf{M^k \ddot u + K^k u = Q^k + S^k} $$

Rather than solving the integrals multiple times in the assembled system, it makes more sense to realize that the integrals themselves do not change. that is, each $M_{ij}^k$ and $K_{ij}^k$ is basically the same. We recall that $N_i^k(x) = \frac{x_{i+1} - x}{h}$ and $N_{i+1}^k(x) = \frac{x-x_i}{h}$. Their derivatives become $N_i'(x) = -\frac{1}{h}$ and $N_{i+1}'(x) = \frac{1}{h}$. In our example $i$ and $j$ are either $i$ or $i+1$, so we find

$$ M_{ii} = \int_{x_i}^{x_{j}} m N_i(x) N_i(x) d\Omega_k = \frac{m}{h^2}\frac{(x_j-x)^3}{3}\vert_{x_i}^{x_{j}} = \frac{mh}{3}$$
$$ M_{ij} = \int_{x_i}^{x_{j}} m N_i(x) N_j(x) d\Omega_k = \frac{m}{h^2}\left[-\frac{(x-x_i)^3}{3}+\frac{h (x-x_i)^2}{2}\right]\vert_{x_i}^{x_{j}} = \frac{mh}{6}$$
$$ M_{ij} = M_{ji} $$
$$ M_{jj} = \int_{x_i}^{x_j} m N_j(x) N_j(x) d\Omega_k = -\frac{m}{h^2}\frac{(x_j-x)^3}{3}\vert_{x_i}^{x_{j}} = \frac{mh}{3}$$

So in total we get
$$ \mathbf{M^k} = \frac{mh}{3} \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix} $$

We do the same for the stiffness matrix as to obtain
$$ \mathbf{K^k} = \frac{EA}{h} \begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix} $$

The internal forces S give:
$$ S_i = 
\sum_{j=i,i+1}^{n_n} \left[ EAN_j'(x) N_i(x) \vert_{x_i}^{x_j}\right] =
EA \left[ N_i'(x)u_i + N_{i+1}'(x) u_{i+1}\right]N_i(x) = 
\frac{EA}{h} (-u_i + u_j) (N_i(x_{i+1})-N_i(x_i))
$$




!!!! CONVERT i and i+1 to a and b


## 4.4 Step 4: Assemble the global system