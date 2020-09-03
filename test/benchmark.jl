using BenchmarkTools
using SPH


function M6_hardcoded(q)
    if q < 0.0
        return 0.0
    elseif q < 1.0
        return (3 - q)^5 - 6*(2 - q)^5 + 15*(1 - q)^5
    elseif q < 2.0
        return (3 - q)^5 - 6*(2 - q)^5
    elseif q < 3.0
        return (3 - q)^5
    else
        return 0.0
    end
end

println("Kernel function benchmarks")
println("Hardcoded M6")
@btime $M6_hardcoded(x) setup=(x=3*rand())

println("Programmatic M6 (and 1st and second derivatives")
@btime $M6[1](x) setup=(x=3*rand())
@btime $(differentiate(M6[1]))(x) setup=(x=3*rand())
@btime $(differentiate(M6[1], 2))(x) setup=(x=3*rand())