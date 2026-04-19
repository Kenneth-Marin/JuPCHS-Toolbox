# Jubuild_PHS
# Julia for port-controlled hamiltonian systems
# By Alejandro Garcés Ruiz
# Reference
# condicion:  todos las matrices continuas en (0,1)
# TO DO: 
# operations that preserves passivity
# |- Paralel 
# |- Feedback
# |- Pre and postmultiplication
# |- Linear transformation
# Plot Hamiltonian in R and R2
# Find the minimum of the Hamiltonian
# Passive PI
# Passive MPC

module JuPCHS

using LinearAlgebra
using DataFrames
using Plots
using ForwardDiff


EPSILON = 1E-6
MAXITER = 500

mutable struct typePCHS    
    J::Function # function J(x,u)
    R::Function # function R(x)
    dH::Function
    H::Function
    G::Function
    nx :: Int64
    nu :: Int64
    dt :: Float64
    xlabel :: Vector
    ulabel :: Vector
    ylabel :: Vector   
    name::String 
end
export BuildPCHS, Simulate, PlotResults, Parallel, Feedback, 
       MatrixTransformation, Heaviside, AdmissibleEquilibrium, 
       Hmin, equilibrium_point, islipschitz

include("JuPCHSConstructors.jl")
include("JuPCHSSimulations.jl")
include("JuPCHSDisplay.jl")
include("JuPCHSOperations.jl") # to do
include("JuPCHSAnalysis.jl") # to do
# JuPCHS_Control.jl
end
