using SPH
using Test, Plots

default(;xlims=(0,3), ylims=(-2, 1), label="")

@testset "SPH.jl" begin
    #@show SPH.Piecewise.p"2x"
    p = plot(M6[1])
    plot!(differentiate(M6[1]))
    plot!(differentiate(M6[1], 2))
    display(p)
    
end
