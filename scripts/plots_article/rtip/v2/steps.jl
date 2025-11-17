include("../../../intro.jl")

f_to = 0.25
polar_amplification = 1.8
f_gmt_pd = 1.2

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
            push!(f[i], params[j + 1, i_f_max] / polar_amplification + f_gmt_pd)
            V[i] = vcat(V[i], ncread(file, "V_sle")[end])
            push!(sr.files[i], file)
        end
    end
end
mutable struct StepRtip{T}
    visc_cases::Vector{String}
    n_visc_cases::Int
    f::Vector{Vector{T}}
    V::Vector{Vector{T}}
    t_end::Vector{Vector{T}}
    files::Vector{Vector{String}}
end
function StepRtip{T}(visc_cases::Vector{String}) where T
    n_visc_cases = length(visc_cases)
    f = [ T[] for _ in 1:n_visc_cases ]
    V = [ T[] for _ in 1:n_visc_cases ]
    t_end = [ T[] for _ in 1:n_visc_cases ]
    files = [ String[] for _ in 1:n_visc_cases ]
    return StepRtip(visc_cases, n_visc_cases, f, V, t_end, files)
end

visc_cases = ["m2stddev", "m1stddev", "nominal", "p1stddev", "p2stddev"]
visc_labels = [L"$-2 \, \sigma$", L"$-1 \, \sigma$", "nominal", L"$+1 \, \sigma$", L"$+2 \, \sigma$"]
sr = StepRtip{Float32}(visc_cases)
prefix = "output/ais/v2/ramps/"
dirs = [
    datadir("$prefix/steps-sigmarange"),
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
        sr.files[i] = sr.files[i][idx]
    end
    return sr
end
sort!(sr)




basins = ["WAIS", "WSB high lat", "RSB", "WSB low lat", "ASB"]
f_gmt_bif = [1.98, 5.87, 7.1, 8.0, 8.65]
V_bif = [54, 42, 31, 20, 13.5]
i_rtip = [[findlast(V .> V_bif[i]) + 1 for i in eachindex(f_gmt_bif)] for V in sr.V]
f_rtip = [sr.f[i][i_rtip[i]] for i in eachindex(i_rtip)]
cls = cgrad(:jet, range(0, stop = 1, length = 6), categorical = true, rev = true)

eqfig = Figure(size=(800, 600), fontsize = 20)
axeq = Axis(eqfig[1, 1])
for i in eachindex(sr.visc_cases)
    lines!(
        axeq,
        sr.f[i],
        sr.V[i],
        label = visc_labels[i];
        linewidth = 4,
        color = cls[i],
    )
    vlines!(axeq, f_rtip[i], color = cls[i], linestyle = :dash)
end
axislegend(axeq, position = :rt)
vlines!(axeq, f_gmt_bif, color = :black, linestyle = :dash, linewidth = 4,
    label = "Bifurcation points")
save(plotsdir("v2/rtip/equil.png"), eqfig)
save(plotsdir("v2/rtip/equil.pdf"), eqfig)

lw = 4
ms = 25
slines_opts = (linewidth = lw, markersize = ms)
set_theme!(theme_latexfonts())
fig = Figure(size=(1000, 500), fontsize = 20)

ax = Axis(fig[1, 1], aspect = AxisAspect(1))
scatterlines!(ax, f_gmt_bif, label = "Bifurcation", color = :black; slines_opts...)
for i in eachindex(basins)
    scatterlines!(ax,
        f_rtip[i],
        label = visc_labels[i],
        color = cls[i];
        slines_opts...,
    )
end
ax.xticks = (eachindex(basins), basins)
ax.xticklabelrotation = pi/8
ax.ylabel = L"$f^\mathrm{crit}$ (K)"
ylims!(ax, (1.5, 9))

ax2 = Axis(fig[1, 2], aspect = AxisAspect(1))
hlines!(ax2, 0, color = :gray70, linestyle = :dash, linewidth = lw)
for i in eachindex(basins)
    scatterlines!(
        ax2,
        f_rtip[i] .- f_gmt_bif,
        label = visc_labels[i],
        color = cls[i];
        slines_opts...,
    )
end
ax2.xticks = (eachindex(basins), basins)
ax2.xticklabelrotation = pi/8
ax2.ylabel = L"$\gamma$ (K)"
ax2.yticks = -1.6:0.2:0.1
ax2.yaxisposition = :right
yl = (-0.8, 0.1)
ylims!(ax2, yl)

Legend(fig[0, 1:2], ax, nbanks = 6)
rowsize!(fig.layout, 0, 5)
colgap!(fig.layout, 5)
fig
save(plotsdir("v2/rtip/gap.png"), fig)
save(plotsdir("v2/rtip/gap.pdf"), fig)