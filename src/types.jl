
"""
    RadialGrid struct for handling radial integration grid
"""
struct RadialGrid
    r::Vector{Float64}    # Radial points
    dr::Vector{Float64}   # Grid spacing
    np::Int               # Number of grid points
end
