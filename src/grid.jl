# module Grid

using LinearAlgebra
using Printf
using Parameters

# export RadialGrid

mutable struct RadialGrid
    zeta::Float64  # nuclear charge
    r::Vector{Float64}  # Radial Points
    rmin::Float64 # Minimum radius
    rmax::Float64   #Maximum radius 
    dr::Vector{Float64} # Grid spacing
    N::Int #Number of points
    # x::Vector{Float64} = Float64[] # x 
    # xmin::Float64 = 1e-4
    # xmax::Float64 = 50.0
    # dx::Vector{Float64} = Float64[]
end

"""
Constructor for RadialGrid
"""
function RadialGrid(;zeta::Float64=1.0, rmin::Float64=1e-4, rmax::Float64=40.0, N::Int=101)

    zeta > 0 || throw(ArgumentError("zeta must be positive"))
    rmin > 0 || throw(ArgumentError("rmin must be positive"))
    rmin < rmax || throw(ArgumentError("rmin must be less than rmax"))
    N >= 101 || throw(ArgumentError("N must be at least 10"))

    xmin, xmax = log(zeta*rmin), log(zeta*rmax)
    x = collect(LinRange(xmin, xmax, N))
    # dx = (x[-1] - x[1]) / N
    dx = x[2] - x[1]
    r = exp.(x)/zeta
    dr = dx * r
    return RadialGrid(zeta, r, rmin, rmax, dr, N)
end

"""
    format_number(x::Number)

Helper function to format numbers in scientific notation when appropriate
"""
function format_number(x::Number)
    if abs(x) < 0.001 || abs(x) > 1000
        return @sprintf("%.3e", x)
    else
        return @sprintf("%.3f", x)
    end
end


function Base.show(io::IO, rg::RadialGrid) 
    println(io,  "RadialGrid:")
    println(io, "  Nuclear charge (ζ) = $(rg.zeta)")
    println(io, "  Grid points (N) = $(rg.N)")
    println(io, "  Range: [$(format_number(rg.rmin)), $(format_number(rg.rmax))] Bohr")
    println(io, "  Grid spacing: $(format_number(rg.dr[1])) to $(format_number(rg.dr[end])) Bohr")
end
#  print(io, rg.zeta, )

# end # Gridd