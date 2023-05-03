# Computational Modelling of Hydraulic and Offshore Structures

## Introduction

The *Computational Modelling* unit is part of all modules B in the *Hydraulic and Offshore Structures* track, embedded in the modules *Marine Renewables (CIEM4220)*, *Dams, Dikes and Breakwaters (CIEM4220)* and *Floating and Submerged Structures (CIEM4230)*. 

The central topic of the unit is the definition and analysis of numerical methods to solve PDEs that govern hydraulic and offshore structures. In this unit it is assumed that the students have already the basic knowledge of the following basic concepts (acquired in the MUDE): 

- Basic numerical methods for ODEs (explicit/implicit). 
- Basic understanding of Finite Differences method for static PDEs. 
- Basic understanding of Finite Element method for static PDEs.

## Unit Learning Objectives

On completion of this unit, the student will be able to
- *LO1:* Construct a conceptual model that represents an hydraulic/offshore engineering application, limited to models that can be constructed from a combination of: point masses, rigid bodies, rods, Euler-Bernoulli beams, geometrically non-linear rods or simple 2-dimensional geometries.
- *LO2*: Apply different numerical methods to solve the equations of motion of the model, subject to typical hydraulic/offshore loads such as: wind, waves and currents.
- *LO3*: Implement the numerical methods and solve the problem using Matlab/Python/others.
- *LO4*: Analyse the results by: validating against analytical solutions or experimental data, identifying the range of applicability of a given method, evaluating errors and assessing the convergence of the solution.

## Forms of instruction

In this unit we will use a problem-based teaching approach. At each week we will provide theoretical material in form of lectures and pre-recorded videos that will later be used in hands-on tutorials and exercises. You are expected to use the concepts practiced during the tutorials and exercises in your final group project. 

Each week we will have a 4 hours workshop where we will cover basic theory (~1h), work on a guided tutorial (~1h) and work on exercises with applications of interest (~2h). In addition, we will have a 2h walk-in feedback session every week. 

**A note on the course structure**

As can be seen here the course is built up of a series of theory pages, as well as tutorial notebooks. It is advised for students to look into the theory pages first. Here a series of texts and videos explains the course material. After that the tutorial pages show ways of implementing the theory lectures. Tutorials are either example-only, with the full example and implementation shown, or they are exercise-solution based. This is indicated in the title. While all exercises have an example solution directly availible at the bottom of the page, it is advised to try it yourself first. In coding there are always different possible paths, so discuss your approach with fellow students or your lecturers. The correctness and efficiency of the solution is key, and even the provided solutions have room for improvement. 

## Assesment on module level

The assessment of this unit will be embedded on the project/exercise portfolio and indivilual exam of each module B. The following percentages will be applied:
-	HOS-B-1 Probabilistic Design: 20%, 
-	HOS-B-2 Numerical Modelling: 25%
-	HOS-B-5 Floating and Submerged Structures: 55%
A successful participation in the module, requires a minimum evaluation of 6.0 for all assessments.

## Lecture and study materials

Lecture slides, lecture notes, papers and chapters from relevant literature, Jupyterbook.

Since this year this Jupyterbook contains the theory, exercises and tutorials from the Computational Modelling unit. For each week you will have the description of the theory, with links to the pre-recorded videos, some solved exercises applying the concepts given in the theory and tutorials with applications to examples relevant to the modules. You will also find a set of exercises to be developed in the workshop sessions. Other than the Python notebooks some comparisons are made using other programs and laguages suchs as [Maple](https://www.maplesoft.com/) and [Julia](https://julialang.org/).  If the installation of additional software is required, we will provide instructions on how to do it.

## Course structure and planning

1. **Week 1**: Review of numerical methods for the solution of Ordinary Differential Equations and Partial Differential Equations.
   1. Taylor series
   2. Forward/Backward Euler
      <br> Look for a [FE](./Module1/w1_t1.ipynb) implementation example.(Module4/w4_t2.ipynb).
2. **Week 2**: Numerical methods for dynamics of systems with multiple rigid bodies. 
   1. Review on derivation of the equation of motion from the Lagrangian (Structural dynamics)
      <br> Look here for a [1DOF](./Module2/w2_t1.ipynb) and [2DOF](./Module2/w2_t2.ipynb) example.
   2. ODE solvers for a system of equations
3. **Week 3**: Numerical methods for dynamic analysis of one-dimensional structures
   1. Dynamics of rods and bars. Finite Differences (FD) and Finite Element Method (FEM) for the transient Laplacian equation
      <br> Look here for an example for a [rod](Module4/w4_t1.ipynb). 
   2. Dynamics of beams. FEM for the transient Euler-Bernoulli beams.
   <br> Look here for an example for a [cantilever beam](Module4/w4_t2.ipynb).
   3. Dynamics of space frame structures. FEM for coupled axial and bending. 
4. **Week 4**: Numerical methods for static analysis of two-dimensional structures
   1. Flow through porous media. FEM for the static Laplacian equation
   2. Settlement of a breakwater. FEM for static linear elasticity
5. **Week 5**: Numerical methods for dynamic analysis of two-dimensional structures
   1. *Dynamic 2D problem*. FEM for dynamic linear elasticity / Stokes flow
6. **Week 6**: Modal analysis
   1. Modal superposition of linear systems
      <br> Look here for a [modal superposition](w5_t1.ipynb) example.
   2. Reduced order models
7. **Week 7**: Numerical methods for geometrically nonlinear structures.
   1. Quasi-static analysis of 1D structures subject to large deformations.
      <br> Look here for [static](w6_t1.ipynb) and [dynamic](w6_t2.ipynb) examples.
   2. Dynamics of 1D structures with large deformations.


## Teaching team ##

::::{grid}
:gutter: 4

:::{grid-item-card} Instructor: J.O. Colomes Gene
```{image} ./images/ProfilePics/Colomes_profilepic.png
:alt: Colo
:class: bg-primary mb-1
:width: 100px
:align: center
```
Mail: [J.O.ColomesGene@tudelft.nl](mailto:J.O.ColomesGene@tudelft.nl)
:::

:::{grid-item-card} Teaching Assistant: J. Modderman
```{image} ./images/ProfilePics/Modderman_profilepic.png
:alt: Modd
:class: bg-primary mb-1
:width: 100px
:align: center
```
Teaching assistant
Mail: [j.modderman@tudelft.nl](mailto:j.modderman@tudelft.nl)
:::

:::{grid-item-card} Teaching Assistant: R. Dekeyser
```{image} ./images/ProfilePics/Dekeyser_profilepic.png
:alt: Deke
:class: bg-primary mb-1
:width: 100px
:align: center
```
Teaching assistant
Mail: [r.dekeyser@student.tudelft.nl](mailto:r.dekeyser@student.tudelft.nl)
:::

::::

## Credits ##

We acknowledge direct contributions to this book from:
- Oriol Colomés
- Jan Modderman
- Ruben Dekeyser

We would also like to akcnowledge contributors from teachers of the course on **Introduction to Computational Dynamics of Offshore Structures** from the master in **Offshore and Dredging Engineering** at TU Delft. The material given in this course has been taken as a basis to build this book.