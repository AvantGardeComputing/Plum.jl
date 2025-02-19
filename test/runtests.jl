using Plum
using Test

@testset "Plum.jl" begin
    # Write your tests here.
end


@testset "grid.jl" begin
    zeta = 1.0
    rmin = 1e-4
    rmax=50.0
    np = 1001
    rg = Plum.RadialGrid(zeta = zeta, rmin = rmin, rmax = rmax, np = np)
    @test rg.dr[1] == 1.3122363377403812e-6
    @test rg.dr[end] == 0.6561181688701899
end
