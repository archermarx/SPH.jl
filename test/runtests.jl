using SPH
using Test

#default(;xlims=(0,3), ylims=(-2, 1), label="")

function M6_hardcoded(q)
    if q < 0.0
        res = 0.0
    elseif q < 1.0
        res = (3 - q)^5 - 6*(2 - q)^5 + 15*(1 - q)^5
    elseif q < 2.0
        res = (3 - q)^5 - 6*(2 - q)^5
    elseif q < 3.0
        res = (3 - q)^5
    else
        res = 0.0
    end
    return res
end

@testset "SPH.jl" begin
   xs = rand(1000) * 6 .- 1
   fs_hardcoded = M6_hardcoded.(xs) .* SPH.σ6[1]
   fs_procedural = M6[1].(xs)
   @test isapprox.(fs_hardcoded, fs_procedural, atol=1e-8) |> all
   #@show xs[a], fs_hardcoded[a], fs_procedural[a]
end
