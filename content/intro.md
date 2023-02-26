# Computational Modelling of Hydraulic and Offshore Structures

## Module Introduction

This Computational modelling unit is part of the *Track module B* in the *Hydraulic and Offshore Structures* track, embedded in the modules *Marine Renewables (CIEM4220)*, *Dams, Dikes and Breakwaters (CIEM4220)* and *Floating and Submerged Structures (CIEM4230)*. 


**Module contents**
The module addresses the design, monitoring and assessment of floating and (semi)-submerged civil engineering solutions for infrastructural and urban development applications, for instance floating bridges and moored tunnels, including control, installation, maintenance and economics. This implies that the structure should be understood as a system as a whole and the varying environmental interactions – aerodynamic, hydrodynamic and soil – be considered in conjunction.
The content of the module contributes to the MUDE and Ethics learning lines.

**Academic and personal skills**
Critical thinking, analytical thinking and problem solving are the main academic skills which will be trained. In addition, students will need to collaborate intensively with their colleagues in a module-wide project so organizational and interpersonal skills will be developed as well.

**Assesment on module level**
The module is assessed at the end of Q4:
-	Group assignment: 60%
-	Individual exam: 40%
Relation between assessment module units (weighing)
-	HOS-B-1 Probabilistic Design: 20%, 
-	HOS-B-2 Numerical Modelling: 25%
-	HOS-B-5 Floating and Submerged Structures: 55%
A successful participation in the module, requires a minimum evaluation of 6.0 for all assessments.

**Prerequisites**
BSc Civil Engineering or Applied Science (or equivalent). MUDE Module, CE Programme Core Module, HOS Track Core Module and HOS Modules A1 or HOS Module A2. 

## Unit introduction

**General information**

**Instructor(s)**
Oriol Colomes, Jeroen Hoving

**Contact Hours / Week (x/x/x/x)**
0/0/0/8

**Participating tracks**
HOS

**Module learning objectives**
MLO-2.
Evaluate hydraulic and offshore structures in interaction with the surrounding environment, using available scientific knowledge.

**Relation to module content**
The central topic of the unit is the definition and analysis of numerical methods to solve PDEs that govern hydraulic and offshore structures. In this unit it is assumed that the students have already the basic knowledge of the following basic concepts (acquired in the MUDE):
-	Basic numerical methods for ODEs (explicit/implicit).
-	Basic understanding of Finite Differences method for static PDEs.
-	Basic understanding of Finite Element method for static PDEs.

**Unit Learning Objectives**
On completion of this unit, the student will be able to
1. ULO-1.	
   Construct a conceptual model that represents an hydraulic/offshore engineering application, limited to models that can be constructed from a combination of: point masses, rigid bodies, rods, Euler-Bernoulli beams, geometrically non-linear rods or simple 2-dimensional geometries.
2. ULO-2.	
   Apply different numerical methods to solve the equations of motion of the model, subject to typical hydraulic/offshore loads such as: wind, waves and currents.
3. ULO-3.	
   Implement the numerical methods and solve the problem using Matlab/Python/others.
4. ULO-4.	
   Analyse the results by: validating against analytical solutions or experimental data, identifying the range of applicability of a given method, evaluating errors and assessing the convergence of the solution.

**Forms of instruction**
Lectures, workshops and design exercise with problem-based teaching. The students are expected to study background material as supplied to them on Brightspace. Assignments are completed outside of lecture and focus on direct application of the theory and methods covered. The design exercise is based on Project-Base-Learning and students will independently develop selected components of common hydraulic and offshore structures. Online and on-campus lectures and on campus workshops.

**Feedback**
Peer feedback, lecturer feedback and self-evaluation.

**Assesment method**
Group assignment, individual exam

**Lecture and study materials**
Lecture slides, lecture notes, papers and chapters from relevant literature, Jupyterbook.

Since this year this Jupyterbook contains the theory, exercises and tutorials from the Computational Modelling unit. For each week you will have the description of the theory, with links to the pre-recorded videos, some solved exercises applying the concepts given in the theory and tutorials with applications to examples relevant to the modules. You will also find a set of exercises to be developed in the workshop sessions. Other than the Python notebooks some comparisons are made using other programs and laguages suchs as [Maple](https://www.maplesoft.com/) and [Julia](https://julialang.org/).  If the installation of additional software is required, we will provide instructions on how to do it.

**Link to the study guide**
(https://studiegids.tudelft.nl/a101_displayCourse.do?course_id=63757)

## Course structure and planning
1. **Week 1**: Review of numerical methods for the solution of Ordinary Differential Equations and Partial Differential Equations.
   1. Forward/Backward Euler
   2. Finite Differences method (FDM)
   3. Finite Elements method (FEM)
2. **Week 2**: Numerical methods for dynamics of systems with multiple rigid bodies. 
   1. Review on derivation of the equation of motion from the Lagrangian (Structural dynamics)
   2. ODE solvers for a system of equations
3. **Week 3**: Numerical methods for dynamic analysis of one-dimensional structures
   1. Dynamics of rods and bars. FEM for the transient Laplacian equation
   2. Dynamics of beams. FEM for the transient Euler-Bernoulli beams.
   3. Dynamics of space frame structures. FEM for coupled axial and bending. 
4. **Week 4**: Numerical methods for static analysis of two-dimensional structures
   1. Flow through porous media. FEM for the static Laplacian equation
   2. Settlement of a breakwater. FEM for static linear elasticity
5. **Week 5**: Numerical methods for dynamic analysis of two-dimensional structures
   1. *Dynamic 2D problem*. FEM for dynamic linear elasticity / Stokes flow
6. **Week 6**: Modal analysis
   1. Modal superposition of linear systems
   2. Reduced order models
7. **Week 6**: Numerical methods for geometrically nonlinear structures.
   1. Quasi-static analysis of 1D structures subject to large deformations.
   2. Dynamics of 1D structures with large deformations.


## Contact details

::::{grid}
:gutter: 4

:::{grid-item-card} Instructor: J.O. Colomes Gene
```{image} ../images/ProfilePics/Colomes_profilepic.png
:alt: Colo
:class: bg-primary mb-1
:width: 100px
:align: center
```
Mail: [J.O.ColomesGene@tudelft.nl](mailto:J.O.ColomesGene@tudelft.nl)
:::

:::{grid-item-card} Teaching Assistant: J. Modderman
```{image} ../images/ProfilePics/Modderman_profilepic.png
:alt: Modd
:class: bg-primary mb-1
:width: 100px
:align: center
```
Teaching assistant
Mail: [j.modderman@tudelft.nl](mailto:j.modderman@tudelft.nl)
:::

:::{grid-item-card} Teaching Assistant: R. Dekeyser
```{image} ../images/ProfilePics/Dekeyser_profilepic.png
:alt: Deke
:class: bg-primary mb-1
:width: 100px
:align: center
```
Teaching assistant
Mail: [r.dekeyser@student.tudelft.nl](mailto:r.dekeyser@student.tudelft.nl)
:::

::::