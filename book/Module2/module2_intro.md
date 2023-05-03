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

# Module 4: 2D Finite Element Modelling

## Introduction

In this first week we will learn how to derive, solve and analyse the equations of motion of rigid body systems.

## Learning Objectives

At the end of this module you will be able to **define and analyze numerical methods to solve the dynamic motion of rigid body systems**. This entails:

1. Characterize a structure as a set of point masses, rigid bodies, rods and beams interacting between each other
2. Define the Equations of Motion of a system through a Hamiltonian approach
3. Define the linearized Equation of Motion of a nonlinear system
4. Define numerical methods to solve a system of ODEs
5. Implement a solver for a system of ODEs
6. Analize and justify the results
