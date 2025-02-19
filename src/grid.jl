# module Grid

using LinearAlgebra
using Printf
# using Parameters

# export RadialGrid

# Add this line at the top to opt out of precompilation
# __precompile__(false)

# @with_kw 
mutable struct RadialGrid
    zeta::Float64  # nuclear charge
    r::Vector{Float64}  # Radial Points
    rmin::Float64 # Minimum radius
    rmax::Float64   #Maximum radius 
    dr::Vector{Float64} # Grid spacing
    np::Int #Number of points
    # x::Vector{Float64} = Float64[] # x 
    # xmin::Float64 = 1e-4
    # xmax::Float64 = 50.0
    # dx::Vector{Float64} = Float64[]
end

"""
Constructor for RadialGrid
"""
function RadialGrid(;
    zeta::Float64=1.0, 
    rmin::Float64=1e-4,
     rmax::Float64=40.0, 
     np::Int=1001)

    zeta > 0 || throw(ArgumentError("zeta must be positive"))
    rmin > 0 || throw(ArgumentError("rmin must be positive"))
    rmin < rmax || throw(ArgumentError("rmin must be less than rmax"))
    np >= 101 || throw(ArgumentError("np must be at least 10"))

    xmin, xmax = log(zeta*rmin), log(zeta*rmax)
    x = collect(LinRange(xmin, xmax, np))
    # dx = (x[-1] - x[1]) / np
    dx = x[2] - x[1]
    r = exp.(x)/zeta
    dr = dx * r
    return RadialGrid(zeta, r, rmin, rmax, dr, np)
end

# use Logarithmic grid from USPP-736 initsubs.f
# a = exp(-aasf) / zed
# b = 1.0 / bbsf
function RadialGrid(
    zeta::Float64=1.0, 
    rmin::Float64=1e-4,
     rmax::Float64=40.0, 
     aasf::Float64=7.0,
     bbsf::Float64=60.0
     )
    zeta > 0 || throw(ArgumentError("zeta must be positive"))
    rmin > 0 || throw(ArgumentError("rmin must be positive"))
    rmin < rmax || throw(ArgumentError("rmin must be less than rmax"))
    aasf > 0 || throw(ArgumentError("aasf must be positive"))
    bbsf > 0 || throw(ArgumentError("bbsf must be positive"))

    # np >= 101 || throw(ArgumentError("np must be at least 10"))

     a = exp(-aasf) / zeta
     b = 1.0 / bbsf

    np::Int = floor(Int, log(rmax / a + 1.0) / b)
    #     ensure that mesh is odd for simpson's rule integration
    np = iseven(np) ? np + 1 : np

    i = collect(LinRange(1, np, np))
    r = a * exp.(b *( i .- 1.0 ) .- 1.0 )
    rab = b* (r .+ a)
    # xmin, xmax = log(zeta*rmin), log(zeta*rmax)
    # x = collect(LinRange(xmin, xmax, np))
    # dx = (x[-1] - x[1]) / np
    # dx = x[2] - x[1]
    # r = exp.(x)/zeta

    dr = rab # dr/di = a * b * e^(b (i-1))
    return RadialGrid(zeta, r, rmin, rmax, dr, np)
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
    println(io, "  Grid points (np) = $(rg.np)")
    println(io, "  Range: [$(format_number(rg.rmin)), $(format_number(rg.rmax))] Bohr")
    println(io, "  Grid spacing: $(format_number(rg.dr[1])) to $(format_number(rg.dr[end])) Bohr")
end


# end # Gridd