include("../../intro.jl")

f_to = 0.25
polar_amplification = 1.8

function aggregate_xp!(sr, dir)
    (;visc_cases, f, V, t_end) = sr
    params = readdlm(joinpath(dir, "info.txt"))
    i_dirs = findfirst(params[1, :] .== "rundir")
    i_f_max = findfirst(params[1, :] .== "hyster.f_max")
    i_visc_case = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
    i_df_dt_max = findfirst(params[1, :] .== "hyster.df_dt_max")
    i_scale = findfirst(params[1, :] .== "isos.viscosity_scaling")
    n_files = size(params, 1) - 1

    for j in 1:n_files
        file = joinpath(dir, "$(params[j + 1, i_dirs])", "yelmo1D.nc")
        visc_case = params[j + 1, i_visc_case]

        if visc_case == "stddev"
            if params[j + 1, i_scale] == -2.0
                visc_case = "m2stddev"
            elseif params[j + 1, i_scale] == -1.0
                visc_case = "m1stddev"
            elseif params[j + 1, i_scale] == 0.0
                visc_case = "nominal"
            elseif params[j + 1, i_scale] == 1.0
                visc_case = "p1stddev"
            elseif params[j + 1, i_scale] == 2.0
                visc_case = "p2stddev"
            end
        end
        i = findfirst(visc_cases .== visc_case)

        t_e = ncread(file, "time")[end]
        if (t_e > 1_000) && (i !== nothing) && (params[j + 1, i_df_dt_max] > 0.5)
            push!(t_end[i], t_e)
            push!(f[i], params[j + 1, i_f_max] / polar_amplification)
            V[i] = vcat(V[i], ncread(file, "V_sle")[end])
        end
    end
end
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

visc_cases = ["m2stddev", "m1stddev", "nominal", "p1stddev", "p2stddev"]
visc_labels = [L"$-2 \, \sigma$", L"$-1 \, \sigma$", "nominal", L"$+1 \, \sigma$", L"$+2 \, \sigma$"]
sr = StepRtip{Float32}(visc_cases)
prefix = "output/ais/ramps/16km"
dirs = [
    datadir("$prefix/steps-sigmarange"),
#     datadir("$prefix/steps-high"),
#     datadir("$prefix/steps-low"),
#     datadir("$prefix/coarse-nomax"),
#     datadir("$prefix/coarse-min"),
]
for dir in dirs
    aggregate_xp!(sr, dir)
end

function Base.sort!(sr::StepRtip{T}) where T
    for i in 1:sr.n_visc_cases
        idx = sortperm(sr.f[i])
        sr.f[i] = sr.f[i][idx]
        sr.V[i] = sr.V[i][idx]
        sr.t_end[i] = sr.t_end[i][idx]
    end
    return sr
end
sort!(sr)

# @show [length(sr.f[i]) for i in 1:sr.n_visc_cases]
# @show [length(sr.V[i]) for i in 1:sr.n_visc_cases]

basins = ["WAIS", "WSB high lat", "RSB", "WSB low lat", "ASB"]
# f_bif = [0.3, 1.8, 2.6, 3.0, 3.5] ./ f_to ./ polar_amplification
# f_bif = [0.4, 2.0, 2.5, 3.0, 3.4] ./ f_to ./ polar_amplification
f_bif = [1.25, 4.3, 6.0, 7.11, 7.89]
V_bif = [55.5, 47.5, 34.5, 23.5, 13.5]
i_rtip = [[findfirst(V .< V_bif[i]) for i in eachindex(f_bif)] for V in sr.V]
f_rtip = [sr.f[i][i_rtip[i]] for i in eachindex(i_rtip)]
# f_bif = f_rtip[1]
@show i_rtip

lw = 4
ms = 25
slines_opts = (linewidth = lw, markersize = ms)
set_theme!(theme_latexfonts())
fig = Figure(size=(1000, 500), fontsize = 20)
cycling_idx = [2, 4, 3, 1, 5]

ax = Axis(fig[1, 1], aspect = AxisAspect(1))
scatterlines!(ax, f_bif, label = "Bifurcation", color = :black; slines_opts...)
for i in eachindex(basins)
    scatterlines!(ax, f_rtip[i], label = visc_labels[i],
        color = Cycled(cycling_idx[i]); slines_opts...)
end
ax.xticks = (eachindex(basins), basins)
ax.xticklabelrotation = pi/8
ax.ylabel = L"$f^\mathrm{crit}$ (K)"
# axislegend(ax, position = :lt, labelsize = 16)
ax.ygridvisible = false
ax.xgridvisible = false
ylims!(ax, (0, 8))

ax_ocean = Axis(fig[1, 1], aspect = AxisAspect(1))
hidexdecorations!(ax_ocean)
ylims!(ax_ocean, (0, 8 * polar_amplification * f_to))
ax_ocean.ylabel = L"$f^\mathrm{crit}_o$ (K)"
ax_ocean.ytickcolor = :gray60
ax_ocean.yticklabelcolor = :gray60
ax_ocean.ylabelcolor = :gray60
ax_ocean.yaxisposition = :right
ax_ocean.ygridvisible = false
ax_ocean.xgridvisible = false

ax2 = Axis(fig[1, 2], aspect = AxisAspect(1))
hlines!(ax2, 0, color = :gray70, linestyle = :dash, linewidth = lw)
for i in eachindex(basins)
    scatterlines!(ax2, f_rtip[i] .- f_bif, label = visc_labels[i],
        color = Cycled(cycling_idx[i]); slines_opts...)
end
ax2.xticks = (eachindex(basins), basins)
ax2.xticklabelrotation = pi/8
ax2.ylabel = L"$\gamma$ (K)"
ax2.ygridvisible = false
ax2.xgridvisible = false
ax2.yticks = -1.6:0.2:0.1
yl = (-1.0, 0.1)
ylims!(ax2, yl)

ax2_ocean = Axis(fig[1, 2], aspect = AxisAspect(1))
hidexdecorations!(ax2_ocean)
ylims!(ax2_ocean, yl .* polar_amplification .* f_to)
ax2_ocean.ylabel = L"$\gamma_\mathrm{o}$ (K)"
ax2_ocean.ytickcolor = :gray60
ax2_ocean.yticklabelcolor = :gray60
ax2_ocean.ylabelcolor = :gray60
ax2_ocean.yaxisposition = :right
ax2_ocean.ygridvisible = false
ax2_ocean.xgridvisible = false

Legend(fig[0, 1:2], ax, nbanks = 6)

rowsize!(fig.layout, 0, 5)
fig
save(plotsdir("16km/rtip/gap.png"), fig)
save(plotsdir("16km/rtip/gap.pdf"), fig)