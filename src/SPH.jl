module SPH

using Piecewise, Plots

export differentiate

include("kernels.jl")

const PARTICLE_MASS = 1.0
const H_INITIAL = 0.1
const ρ_INITIAL = 0.5
const T_INITIAL = 1

# TODO: Generalize for n dimensions
struct Particle{S}
    x::S
    y::S
    ρ::S
    h::S
    vx::S
    vy::S
    T::S
    Particle(x::S, y::S, ρ::S, h::S, vx::S, vy::S, T::S) where S = new{S}(x, y, ρ, h, vx, vy, T)
end

Particle(x::T, y::T) where T = Particle(x, y, ρ_INITIAL * one(T), H_INITIAL)

Particle(x::T, y::T, ρ::T, h::T) where T = Particle(x, y, ρ, h, zero(T), zero(T), one(T))


particles_in_box(n, w=1, h=1) = [Particle(w*rand(), h*rand()) for i in 1:n]

function plot_particles(particles; kwargs...)
    xs = [p.x for p in particles]
    ys = [p.y for p in particles]
    hs = [p.h for p in particles]
    ρs = [p.ρ for p in particles]

    p = plot()
    #draw_circles!(xs, ys, hs, lw = 0.5, fillalpha = 0)
    scatter!(xs, ys; mz = ρs, markerstrokewidth=0, markersize = 20, kwargs...)
    #surface!(xs, ys, ρs, camera=(0, 90); kwargs...)
    return p
end

function draw_circles!(xs, ys, hs; kwargs...)
    θ = LinRange(0, 2π, 50)
    for (x, y, h) in zip(xs, ys, hs)
        circle_x = x .+ h*sin.(θ)
        circle_y = y .+ h*cos.(θ)
        plot!(circle_x, circle_y; seriestype = [:shape], kwargs...)
    end
end

const MAX_TIME_STEPS = 100

# TODO: Use sparse upper triangular matrix for this
compute_distances(particles) = [hypot(p.x - q.x, p.y - q.y) for p in particles, q in particles]

function update_variable(particles, kernel = M6; binary = false, η = 1, dimension = 2)
    distances = compute_distances(particles)
    n_particles = length(particles)

    #new_particles = Vector{Particle}(undef, n_particles)
    new_particles = copy(particles)

    w  = kernel[dimension]
    ∂w = differentiate(kernel[dimension])

    d = dimension
    m = PARTICLE_MASS

    #binary = true
    
    for i in 1:n_particles
        xᵢ, yᵢ  = particles[i].x, particles[i].y
        ρᵢ, hᵢ  = particles[i].ρ, particles[i].h

        @views rs = distances[:, i]
        res = 1
        niter = 1
       
        if binary
            maxh = 1000
            minh = 0
            hᵢ = (maxh + minh)/2
        end

        q = zeros(n_particles)
        Wᵢ = zeros(n_particles)
        ∂Wᵢ = zeros(n_particles)

        #while abs(res) > 1e-8 && niter < MAX_TIME_STEPS
        #    q = rs / hᵢ
        #    Wᵢ = w.(q) / hᵢ^d
        #    ∂Wᵢ = ∂w.(q) / hᵢ^d
        #
        #    f  = m*((η/hᵢ)^d - sum(Wᵢ))
        #
        #    if binary
        #        if signbit(f)
        #            maxh = hᵢ
        #        else
        #            minh = hᵢ
        #        end
        #        newh = (maxh + minh)/2
        #        res = newh - hᵢ
        #        hᵢ = newh
        #    else
        #        f′ = -m*(d * η^d / hᵢ^(1 + d) + sum(q/hᵢ .* ∂Wᵢ .+ d/hᵢ * Wᵢ))    
        #        res = - f/f′
        #        hᵢ += res
        #    end
        #end

        #ρᵢ = m*(η/hᵢ)^d
        qs = rs / particles[i].h
        ρᵢ = sum(m) 

        Wᵢ = w.(q) / hᵢ^d
        #    ∂Wᵢ = ∂w.(q) / hᵢ^d

        #dhdρ = -hᵢ/ρᵢ/d

        dWdh = q/h * ∂Wᵢ - d/h * Wᵢ

        Ωᵢ = 1 - dhdρ * m * sum(dWdh)

        Pᵢ = ρᵢ * particles[i].T

        dvdt = 0

        #@show hᵢ, ρᵢ

        new_particles[i] = Particle(xᵢ, yᵢ, ρᵢ, hᵢ, vx, vy, T)
    end
    return new_particles
end

function update(particles; dt = 0.01, dimension = 2, kernel = M6, η = 1)
    distances = compute_distances(particles)
    n_particles = length(particles)

    #new_particles = Vector{Particle}(undef, n_particles)
    new_particles = copy(particles)

    w  = kernel[dimension]
    ∂w = differentiate(kernel[dimension])

    d = dimension
    m = PARTICLE_MASS

    h = 0.1
    q = distances / h

    # Update density
    new_ρ = zeros(n_particles)

    for i = 1:n_particles
        ρᵢ = 0.0
        for j in 1:n_particles
            if i != j
                Wᵢⱼ = 1/h^d * w(q[i, j])
                ρᵢ += Wᵢⱼ
            end
        end

        new_ρ[i] = ρᵢ
    end

    for i in 1:n_particles
        ρᵢ = 0.0
        dvx_dt = 0.0
        dvy_dt = 0.0

        particle_i = particles[i]
        xᵢ = particle_i.x
        yᵢ = particle_i.y

        for j in 1:n_particles
            if i != j
                ∂Wᵢⱼ = 1/h^d * ∂w(q[i, j])

                particle_j = particles[j]
                xⱼ = particle_j.x
                yⱼ = particle_j.y

                ρᵢ = new_ρ[i]
                ρⱼ = new_ρ[j]

                pᵢ = ρᵢ * particle_i.T
                pⱼ = ρⱼ * particle_j.T

                term = ∂Wᵢⱼ/distances[i, j]/h

                ∇Wx = term * (xᵢ - xⱼ)
                ∇Wy = term * (yᵢ - yⱼ)

                term = (pᵢ/ρᵢ^2 + pⱼ/ρⱼ^2)

                dvx_dt -= term * ∇Wx
                dvy_dt -= term * ∇Wy
                
                #isnan(dvx_dt) && println("dvx_dt nan at $i, $j")
                #isnan(dvy_dt) && println("dvy_dt nan at $i, $j")
            end
        end
        ρᵢ = m*ρᵢ
        dvx_dt = m*dvx_dt
        dvy_dt = m*dvy_dt
        vx = particle_i.vx + dt * dvx_dt
        vy = particle_i.vy + dt * dvy_dt
        
        x = particle_i.x + vx*dt
        y = particle_i.y + vy*dt

        if x >= 1
            x -= 1
        elseif x <= 0
            x += 1
        elseif y > 1
            y -= 1
        elseif y <= 0
            y += 1
        end

        new_particles[i] = Particle(x, y, ρᵢ, h, vx, vy, particle_i.T)
    end

    return new_particles
end

end