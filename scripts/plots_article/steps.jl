include("../intro.jl")

mutable struct StepRtip{T}
    visc_cases::Vector{String}
    n_visc_cases::Int
    f::Vector{Vector{T}}
    V::Vector{Vector{T}}
    t_end::Vector{Vector{T}}
end
function StepRtip{T}(visc_cases::Vector{String}) where T
    n_visc_cases = length(visc_cases)
    f = [ T[] for _ in 1:n_visc_cases ]
    V = [ T[] for _ in 1:n_visc_cases ]
    t_end = [ T[] for _ in 1:n_visc_cases ]
    return StepRtip(visc_cases, n_visc_cases, f, V, t_end)
end

visc_cases = ["min", "-2std", "scale", "+2std", "max"]
visc_labels = ["MIN", "LOW", "SCALE", "HIGH", "MAX"]
sr = StepRtip{Float32}(visc_cases)
dirs = [
    datadir("output/ais/ramps/16km/steps-minomax"),
    datadir("output/ais/ramps/16km/steps-min"),
]
function aggregate_xp!(sr, dir)
    (;visc_cases, f, V, t_end) = sr
    params = readdlm(joinpath(dir, "info.txt"))
    i_dirs = findfirst(params[1, :] .== "rundir")
    i_f_max = findfirst(params[1, :] .== "hyster.f_max")
    i_visc_case = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
    n_files = size(params, 1) - 1
    f_to = 0.25

    for j in 1:n_files
        file = joinpath(dir, "$(params[j + 1, i_dirs])", "yelmo1D.nc")
        visc_case = params[j + 1, i_visc_case]
        i = findfirst(visc_cases .== visc_case)
        t_e = ncread(file, "time")[end]
        if t_e > 10_000
            push!(t_end[i], t_e)
            push!(f[i], params[j + 1, i_f_max] * f_to)
            V[i] = vcat(V[i], ncread(file, "V_sle")[end])
        end
    end
end

for dir in dirs
    aggregate_xp!(sr, dir)
end

# @show [length(sr.f[i]) for i in 1:sr.n_visc_cases]
# @show [length(sr.V[i]) for i in 1:sr.n_visc_cases]

basins = ["Aurora", "Recovery", "Aurora+Recovery", "Wilkes"]
f_bif = [1.8, 2.6, 3.0, 3.5]
V_bif = [46, 35, 27, 14.5]
i_rtip = [[findfirst(V .< V_bif[i]) for i in eachindex(f_bif)] for V in sr.V]
@show i_rtip
f_rtip_3 = sr.f[3][i_rtip[3]]
f_rtip_5 = sr.f[5][i_rtip[5]]

lw = 4
ms = 25
slines_opts = (linewidth = lw, markersize = ms)
set_theme!(theme_latexfonts())
fig = Figure(size=(500, 500), fontsize = 20)
ax = Axis(fig[1, 1], aspect = AxisAspect(1))
scatterlines!(ax, f_bif, label = "Bifurcation"; slines_opts...)
scatterlines!(ax, f_rtip_3, label = "NOM"; slines_opts...)
scatterlines!(ax, f_rtip_5, label = "MAX"; slines_opts...)
ax.xticks = (eachindex(basins), basins)
ax.xticklabelrotation = pi/8
ax.ylabel = "Critical ocean forcing (K)"
axislegend(ax, position = :lt, labelsize = 16)
fig
