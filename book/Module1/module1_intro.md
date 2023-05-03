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

# Module 1: Numerical methods for ODEs

## Introduction 

In this chapter we will go over the contents given in the *Modelling, Uncertainty and Data for Engineers (MUDE)* module that are relevant for this unit. These are the following:
1. Taylor Series Expansions (week 1.5)
2. Solution to Ordinary Differential Equations (week 1.5)

## Learning Objectives

At the end of this module you will be able to:
1. Define and analyse numerical methods to solve ODEs.
  1. Define a simple solver to approximate solutions of ODEs based on Taylor Series
  2. Quantify the numerical error of an approximated solution
  3. Define adaptive time stepping approaches to control the numerical error
  4. Distinguish between different ODE solvers