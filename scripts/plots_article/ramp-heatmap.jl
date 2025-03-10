include("../intro.jl")

mutable struct HeatmapRtip{T}
    visc_cases::Vector{String}
    n_visc_cases::Int
    f::Vector{Vector{T}}
    dfdt::Vector{Vector{T}}
    V::Vector{Vector{T}}
    t_end::Vector{Vector{T}}
    paths::Vector{Vector{String}}
end
function HeatmapRtip{T}(visc_cases::Vector{String}) where T
    n_visc_cases = length(visc_cases)
    f = [ T[] for _ in 1:n_visc_cases ]
    dfdt = [ T[] for _ in 1:n_visc_cases ]
    V = [ T[] for _ in 1:n_visc_cases ]
    t_end = [ T[] for _ in 1:n_visc_cases ]
    paths = [ String[] for _ in 1:n_visc_cases ]
    return HeatmapRtip(visc_cases, n_visc_cases, f, dfdt, V, t_end, paths)
end

function aggregate_xp!(hr::HeatmapRtip, dir::String)

    (;visc_cases, f, dfdt, V, t_end, paths) = hr
    params = readdlm(joinpath(dir, "info.txt"))
    i_dirs = findfirst(params[1, :] .== "runid")
    i_df_dt_max = findfirst(params[1, :] .== "hyster.df_dt_max")
    i_f_max = findfirst(params[1, :] .== "hyster.f_max")
    i_visc_case = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
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
        if (t_e > 10 && i !== nothing)
            push!(t_end[i], t_e)
            push!(f[i], params[j + 1, i_f_max])
            push!(dfdt[i], params[j + 1, i_df_dt_max])
            push!(V[i], ncread(file, "V_sle")[end])
            push!(paths[i], joinpath(dir, "$(params[j + 1, i_dirs])"))
        end
    end
end

visc_cases = ["m2stddev", "m1stddev", "nominal", "p1stddev", "p2stddev"]
visc_labels = [L"$-2 \, \sigma$", L"$-1 \, \sigma$", "nominal", L"$+1 \, \sigma$", L"$+2 \, \sigma$"]
hr = HeatmapRtip{Float32}(visc_cases)

region = "recovery"
if region == "recovery"
    dirs = [
        datadir("output/ais/ramps/16km/ramps-recovery-sigmarange"),
    ]
else
    dirs = [
        datadir("output/ais/ramps/16km/ramps-sigmarange"),
    ]
end

for dir in dirs
    aggregate_xp!(hr, dir)
end
@show [length(f) for f in hr.f]


# Compute SSP slopes
polar_amplification = 1.8
ssps = load_ssps(polar_amplification; wrt = :pd)
mean_slope(t, x) = [ (x[i] - x[1])/(t[i] - t[1]) for i in 2:length(x) ]
mean_slope(x) = mean_slope(view(x, :, 1), view(x, :, 2))
dssps = [mean_slope(ssp) for ssp in ssps]
ssp_labels = ["SSP1", "SSP2", "SSP3", "SSP5"]
ssp_colors = [:darkblue, :lightblue, :orange, :darkred]
dfdt_min = minimum([minimum(dssp) for dssp in dssps])
dfdt_max = maximum([maximum(dssp) for dssp in dssps])

# Compute values related to 2D maps
k = 3

if region == "recovery"
    i1, j1 = 110, 190
    d1, d2 = 130, 115
    f_hm = [10.8, 11.0]
    dfdt_hm = [1e-2, 3e-2]
    f_hm_jet = collect(10.2:0.2:11.4)
    dfdt_hm_jet = [3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]
    f_bif = 10.7
    cmap = (colormap = cgrad([:lightcoral, :white, :cornflowerblue]), colorrange = (31, 36),
        lowclip = :lightcoral, highclip = :cornflowerblue)
else
    i1, j1 = 60, 105
    d1, d2 = 130, 115
    f_hm = [2.0, 2.2]
    dfdt_hm = [1e-2, 3e-2]
    f_hm_jet = collect(1.4:0.2:2.6)
    dfdt_hm_jet = [3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]
    f_bif = 2.1
    cmap = (colormap = cgrad([:lightcoral, :white, :cornflowerblue]), colorrange = (53, 57),
        lowclip = :lightcoral, highclip = :cornflowerblue)
end

idx = [argmin( abs.(hr.f[k] .- f_hm[i]) + abs.(hr.dfdt[k] .- dfdt_hm[i]) ) for i in 1:2]
idx_jet = [argmin( abs.(hr.f[k] .- f_hm_jet[i]) + abs.(hr.dfdt[k] .- dfdt_hm_jet[i]) )
    for i in eachindex(f_hm_jet)]

set_theme!(theme_latexfonts())
fs = 24
ms1 = 5
ms2 = 20
axasp = 1.2
lw = 3
rw_cbar = 0.6
hlmarker = '♥'
fig = Figure(size = (1100, 800), fontsize = fs)

jetmap = cgrad(:copper, range(0, stop = 1, length = 8), categorical = true, rev = true)
framecolors = [:purple4, :purple2]

ax_ramp = Axis(fig[2, 1], aspect = AxisAspect(axasp), valign = :top)
ax_ramp.ylabel = "Forcing (K)"
ax_ramp.xlabel = "Time (kyr)"
ax_ramp.xaxisposition = :top
ax_ramp.xticks = 1:2
ax_ramp.xminorticks = 0.5:0.5:2.5
ax_ramp.xminorticksvisible = true
ax_ramp.xminorgridvisible = true
xlims!(ax_ramp, (0.92, 2.1))

ax_vol = Axis(fig[3, 1], aspect = AxisAspect(axasp), valign = :top)
ax_vol.ylabel = L"$V_\mathrm{AIS}$ (m SLE)"
ax_vol.xticklabelsvisible = true
ax_vol.xticksvisible = true
# ax_vol.xlabel = "Time (kyr)"
ax_vol.xticks = 0:5:15
ax_vol.yticks = 0:10:60
ax_vol.xminorticks = 0:2.5:20
ax_vol.xminorticksvisible = true
xlims!(ax_vol, (0, 20))

ax_hm = [Axis(fig[2:3, 2:3], aspect = DataAspect()), Axis(fig[2:3, 4:5], aspect = DataAspect())]
axs = [Axis(fig[4, j], aspect = AxisAspect(1)) for j in 1:hr.n_visc_cases]


logxticks = (-4:2:0, [L"10^{%$n} $\,$" for n in -4:2:0])
axs[1].ylabel = "Max forcing (K)"
# axs[hr.n_visc_cases].ylabel = "Max forcing (K)"
axs[hr.n_visc_cases].yaxisposition = :right

for j in eachindex(hr.visc_cases)
    
    if length(hr.f[j]) > 1
        heatmap!(axs[j], log10.(hr.dfdt[j]), hr.f[j], hr.V[j]; cmap...)
        scatter!(axs[j], log10.(hr.dfdt[j]), hr.f[j], markersize = ms1, color = :white)
        hlines!(axs[j], [f_bif], color = :gray10, linestyle = :dash, linewidth = 4)
    end
    vlines!(axs[j], log10.([dfdt_min, dfdt_max]), color = :darkred, linewidth = 3)

    axs[j].title = visc_labels[j]
    # j in (1,3,5) ? axs[j].xlabel = L"$\mathrm{log}_{10}$ slope (K)" : nothing
    axs[j].xlabel = L"$\mathrm{log}_{10}$ slope (1)"
    axs[j].xticks = -4:1:0
    axs[j].xminorticks = -4:0.5:0
    axs[j].xminorticksvisible = true
    if region == "recovery"
        axs[j].yticks = 10.0:0.4:11.4
    else
        axs[j].yticks = 1.4:0.4:2.6
    end

    if region == "recovery"
        ylims!(axs[j], (9.9, 11.5))
        xlims!(axs[j], (log10(0.0003)-0.25, log10(0.3)+0.25))
    else
        ylims!(axs[j], (1.3, 2.7))
        xlims!(axs[j], (-4.2, -0.3))
    end

    if 1 < j # < hr.n_visc_cases
        axs[j].yticksvisible = false
        axs[j].yticklabelsvisible = false
    end
end


for i in eachindex(f_hm_jet)
    scatter!(axs[k], log10.(hr.dfdt[k][idx_jet[i]]), hr.f[k][idx_jet[i]], markersize = ms2,
        color = jetmap[i], marker = hlmarker)
    file1D = joinpath(hr.paths[k][idx_jet[i]], "yelmo1D.nc")
    lines!(ax_ramp, ncread(file1D, "time") ./ 1e3, ncread(file1D, "hyst_f_now"),
        color = jetmap[i], linewidth = 3)
    lines!(ax_vol, ncread(file1D, "time") ./ 1e3, ncread(file1D, "V_sle"),
        color = jetmap[i], linewidth = 3)
end

# The inset axis
X = ncread(joinpath(hr.paths[k][1], "yelmo2D.nc"), "x2D")
nx, ny = size(X)
XX, YY = ndgrid(1:nx, 1:ny)
dcrop = 10
ii = dcrop+1:nx-dcrop
jj = dcrop+1:ny-dcrop
inset_axs = [
    Axis(fig[2:3, i],
    width=Relative(0.3),
    height=Relative(0.3),
    halign=0.97,
    valign=0.11,
    backgroundcolor=:white) for i in [2:3, 4:5]]
[hidedecorations!(inset_ax) for inset_ax in inset_axs]


for i in eachindex(f_hm)
    hidedecorations!(ax_hm[i])
    file = joinpath(hr.paths[k][idx[i]], "yelmo2D.nc")
    time = ncread(file, "time")
    nt = length(time)
    H_ice = ncread(file, "H_ice", start = [i1, j1, nt], count = [d1, d2, 1])[:, :, 1]
    z_bed = ncread(file, "z_bed", start = [i1, j1, nt], count = [d1, d2, 1])[:, :, 1]
    heatmap!(ax_hm[i], z_bed; cmaps["z_bed"]...)
    heatmap!(ax_hm[i], H_ice + z_bed; cmaps["z_srf"]...)
    scatter!(axs[k], log10.(hr.dfdt[k][idx[i]]), hr.f[k][idx[i]], markersize = ms2,
        color = framecolors[i], marker = hlmarker)
    file1D = joinpath(hr.paths[k][idx[i]], "yelmo1D.nc")
    contour!(ax_hm[i], ncread(file, "f_grnd", start = [i1, j1, 1], count = [d1, d2, 1])[:, :, 1] .> 1e-8,
        color = :black, linewidth = 2)
    contour!(ax_hm[i], ncread(file, "f_grnd", start = [i1, j1, nt], count = [d1, d2, 1])[:, :, 1] .> 1e-8,
    color = :red, linewidth = 2)
    lines!(ax_ramp, ncread(file1D, "time") ./ 1e3, ncread(file1D, "hyst_f_now"),
        color = framecolors[i], linewidth = lw)
    lines!(ax_vol, ncread(file1D, "time") ./ 1e3, ncread(file1D, "V_sle"),
        color = framecolors[i], linewidth = lw)
    heatmap!(inset_axs[i], ncread(joinpath(hr.paths[k][1], "yelmo2D.nc"), "H_ice", start = [1, 1, 1], count = [-1, -1, 1])[ii, jj, 1] .> 1e-8, colorrange = (1e-8, 1), colormap = cgrad([:white, :gray70]), lowclip = :white, highclip = :gray70)
    contour!(inset_axs[i], ((i1 .<= XX .<= i1+d1) .&& (j1 .<= YY .<= j1+d2))[ii, jj],
        color = :black, linewidth = 2)
    text!(inset_axs[i], 10, 10, font = :bold, color = :black, fontsize = fs-4,
        text = "t = $(Int(round(time[nt] / 1e3, digits = 0))) kyr")
    ax_hm[i].leftspinecolor = framecolors[i]
    ax_hm[i].rightspinecolor = framecolors[i]
    ax_hm[i].topspinecolor = framecolors[i]
    ax_hm[i].bottomspinecolor = framecolors[i]
    ax_hm[i].spinewidth = 5
end



Colorbar(fig[1, 2:3], label = "Bed elevation (km)", vertical = false,
    width = Relative(rw_cbar), flipaxis = true, ticks = latexifyticks(-6:2:4, 1e3),
    halign = :left; cmaps["z_bed"]...)
Colorbar(fig[1, 3:4], label = "Ice surface elevation (km)", vertical = false,
    width = Relative(rw_cbar), flipaxis = true, ticks = (vcat(1, 1e3:1e3:4e3),
    string.(0:4)); cmaps["z_srf"]...)
Colorbar(fig[1, 4:5], label = L"$V_\mathrm{AIS}$ (m SLE)", vertical = false,
    width = Relative(rw_cbar), flipaxis = true, halign = :right; cmap...)


colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)
# rowgap!(fig.layout, 1, -10)
# rowgap!(fig.layout, 2, -10)
# rowgap!(fig.layout, 3, -40)

rowgap!(fig.layout, 1, -60)
rowgap!(fig.layout, 2, -20)
rowgap!(fig.layout, 3, -20)
# rowgap!(fig.layout, 3, -30)
rowsize!(fig.layout, 1, 10)

fig
save(plotsdir("16km/rtip/ramp-heatmap-$region.png"), fig)