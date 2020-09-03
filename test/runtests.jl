using SPH
using Test, Plots

default(;color = :thermal, aspect_ratio = 1, label = "", xaxis=("x", (-0.1, 1.1)), yaxis = ("y", (-0.1, 1.1)))

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
   @test isapprox.(fs_hardcoded, fs_procedural, atol=1e-12) |> all
   #@show xs[a], fs_hardcoded[a], fs_procedural[a]
end

@testset "Density calculations" begin

    particles = SPH.particles_in_box(1)
    @test SPH.compute_distances(particles) ≈ [0.0]

    particles = [SPH.Particle(0.0, 0.0), SPH.Particle(0.0, 1.0)]
    @test SPH.compute_distances(particles) ≈ [0.0 1.0; 1.0 0.0]

    push!(particles, SPH.Particle(1.0, 1.0))
    @test SPH.compute_distances(particles) ≈ [0.0 1.0 √2; 1.0 0.0 1.0; √2 1.0 0.0]

end

N = 100
particles = SPH.particles_in_box(N)

NEdge = floor(Int, sqrt(N))
xs = LinRange(0, 1, NEdge)
ys = LinRange(0, 1, NEdge)

#left = SPH.Particle.(zeros(NEdge), ys)
#bottom = SPH.Particle.(xs, zeros(NEdge))
#top = SPH.Particle.(xs, ones(NEdge))
#right = SPH.Particle.(ones(NEdge), ys)

#append!(particles, left)
#append!(particles, right)
#append!(particles, top)
#append!(particles, bottom)


timestep = 0.01
anim = @animate for i in 1:500
    global particles = SPH.update(particles, dt=timestep);
    if i % 2 == 0
        SPH.plot_particles(particles, clims = (0, 120))
    end
end
gif(anim, "test.gif", fps = 25)

