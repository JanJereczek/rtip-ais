
@kwdef struct PowerLaw{T<:AbstractFloat}
    q::T = 0.75
    u0::T = 100
end

@kwdef struct RegularizedCoulombLaw{T<:AbstractFloat}
    q::T = 0.2
    u0::T = 100
end

friction(u::T, law::PowerLaw{T}) where T = u / ( law.u0^law.q + u^(1-law.q) )
friction(u::T, law::RegularizedCoulombLaw{T}) where T = (u / (law.u0 + u))^law.q

u_vec = collect(range(0, stop = 100, step = 0.1))
T = eltype(u_vec)
pl = PowerLaw{T}()
rcl = RegularizedCoulombLaw{T}()

f_power = map(u -> friction(u, pl), u_vec)
f_coulomb = map(u -> friction(u, rcl), u_vec)

set_theme!(theme_latexfonts())
fig = Figure(fontsize = 20)
ax = Axis(fig[1, 1], xlabel = "Sliding velocity (m/yr)", ylabel = "Friction coefficient (1)")
lines!(ax, u_vec, f_power, label = "Power law", linewidth = 3)
lines!(ax, u_vec, f_coulomb, label = "Regularized Coulomb law", linewidth = 3)
save(plotsdir("friction_laws.png"), fig)