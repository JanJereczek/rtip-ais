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
visc_num_labels = [L"$\textbf{e} \quad -2 \, \sigma$", L"$\textbf{f} \quad -1 \, \sigma$", L"\textbf{g} \quad nominal $\,$",
    L"$\textbf{h} \quad +1 \, \sigma$", L"$\textbf{i} \quad +2 \, \sigma$"]
hr = HeatmapRtip{Float32}(visc_cases)

region = "recovery"
if region == "recovery"
    dirs = [
        datadir("output/ais/ramps/16km/ramps-recovery-sigmarange"),
        datadir("output/ais/ramps/16km/ramps-recovery-sigmarange-extlow"),
        datadir("output/ais/ramps/16km/ramps-recovery-sigmarange-long"),
    ]
elseif region == "wais"
    dirs = [
        datadir("output/ais/ramps/16km/ramps-sigmarange"),
    ]
end

for dir in dirs
    aggregate_xp!(hr, dir)
end
@show [length(f) for f in hr.f]


# Compute SSP slopes
f_pd = 1.2
f_to = 0.25
polar_amplification = 1.8
pa = polar_amplification
ssps = load_ssps(polar_amplification; wrt = :pd)
mean_slope(t, x) = [ (x[i] - x[1])/(t[i] - t[1]) for i in 2:length(x) ]
mean_slope(x) = mean_slope(view(x, :, 1), view(x, :, 2))
dssps = [mean_slope(ssp) for ssp in ssps]
ssp_labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
ssp_colors = [:darkblue, :lightblue, :darkorange, :darkred]
dfdt_min = minimum([minimum(dssp) for dssp in dssps])
dfdt_max = maximum([maximum(dssp) for dssp in dssps])
dfdt_min = 10 ^ floor(log10(dfdt_min))
dfdt_max = 10 ^ ceil(log10(dfdt_max)) # Let's take some margin

# Compute values related to 2D maps
k = 3
if region == "recovery"
    i1, i2 = 110, 190
    d1, d2 = 130, 115
    f_hm = [10.2, 10.4]
    dfdt_hm = [1e-3, 3e-3]
    f_hm_jet = collect(10.0:0.2:11.2)
    dfdt_hm_jet = [3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1]
    f_bif = 10.7
    cmap = (colormap = cgrad([:lightcoral, :white, :cornflowerblue]), colorrange = (31, 36),
        lowclip = :lightcoral, highclip = :cornflowerblue)
elseif region == "wais"
    i1, i2 = 60, 105
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
fs = 20
ms1 = 5
ms2 = 15
axasp = 1.15
lw = 3
rw_cbar = 0.6
viscmarker = :star5
fig = Figure(size = (1100, 770), fontsize = fs)

jetmap = cgrad(:copper, range(0, stop = 1, length = 8), categorical = true, rev = true)
framecolors = [:seagreen, :purple2]

ra = fig[2, 1] = GridLayout()
ha = fig[3, 1:3] = GridLayout()
ax_ramp1 = Axis(ra[1, 1], valign = :top)#, aspect = AxisAspect(axasp)
ax_ramp2 = Axis(ra[1, 2], valign = :top)
ax_ramp1.ylabel = L"GMT anomaly, $f$ (K)"
ax_ramp1.xlabel = "           Time (kyr)"
ax_ramp1.xaxisposition = :top
ax_ramp2.xaxisposition = :top
# ax_ramp1.xminorticksvisible = true
# ax_ramp1.xminorgridvisible = true

ax_vol = Axis(ra[2, 1:2], valign = :top)#, aspect = AxisAspect(axasp)
ax_vol.ylabel = L"$V_\mathrm{AIS}$ (m SLE)"
ax_vol.xticklabelsvisible = true
ax_vol.xticksvisible = true

ax_ramp1.xticks = 1:2
ax_ramp2.xticks = 10:10:45
ax_ramp1.xminorticks = 0.5:0.5:2.5
ax_ramp1.xminorticks = 5:5:45
xlims!(ax_ramp1, (0.92, 2.1))
xlims!(ax_ramp2, (2.1, 45))
ax_vol.xminorticks = 0:5:60
ax_ramp2.yticklabelsvisible = false
ax_ramp2.yticksvisible = false

if region == "recovery"
    ylims!(ax_ramp1, (-0.1, 6.8) .+ f_pd)
    ylims!(ax_ramp2, (-0.1, 6.8) .+ f_pd)
    ax_ramp1.yticks = 1:2:7
    ax_vol.xticks = 0:20:60
    ax_vol.xminorticksvisible = true
    xlims!(ax_vol, (0, 45))
elseif region == "wais"
    ylims!(ax_ramp1, (-0.1, 1.6) .+ f_pd)
    ylims!(ax_ramp2, (-0.1, 1.6) .+ f_pd)
    ax_ramp1.yticks = 1:0.5:3
    ax_vol.xticks = 0:10:30
    ax_vol.yticks = 0:2:60
    xlims!(ax_vol, (0, 30))
end

ax_hm = [Axis(fig[2, 2], aspect = DataAspect()), Axis(fig[2, 3], aspect = DataAspect())]
axs = [Axis(ha[1, j], aspect = AxisAspect(1)) for j in 1:hr.n_visc_cases]

logxticks = (-4:2:0, [L"10^{%$n} $\,$" for n in -4:2:0])
axs[1].ylabel = L"Max GMT anomaly, $f^\mathrm{max}$ (K)"
axs[hr.n_visc_cases].yaxisposition = :right
rsl_index = argmin( (hr.f[1] .- 10.2).^2 .+ (hr.dfdt[1] .- 3e-2).^2 )
file_rsl_ref = joinpath(hr.paths[1][rsl_index], "yelmo2D.nc")
X, Y = ncread(file_rsl_ref, "x2D"), ncread(file_rsl_ref, "y2D")
x, y = ncread(file_rsl_ref, "xc"), ncread(file_rsl_ref, "yc")
x1, x2, y1, y2 = -500, -100, 900, 1400
bbox = (x1 .< X .< x2) .&& (y1 .< Y .< y2)

visc_colors = cgrad(:jet, range(0, stop = 1, length = hr.n_visc_cases + 1), categorical = true,
    rev = true)
visc_colors = [visc_colors[1], visc_colors[2], :darkgreen, visc_colors[4], visc_colors[5]]

for j in eachindex(hr.visc_cases)
    
    if length(hr.f[j]) > 1
        heatmap!(axs[j], log10.(hr.dfdt[j] ./ pa), hr.f[j] ./ pa .+ f_pd, hr.V[j]; cmap...)
        scatter!(axs[j], log10.(hr.dfdt[j] ./ pa), hr.f[j] ./ pa .+ f_pd, markersize = ms1,
            color = :white)
        if region == "recovery"
            scatter!(axs[j], log10.(hr.dfdt[j][rsl_index] ./ pa),
                hr.f[j][rsl_index] ./ pa .+ f_pd,
                markersize = ms2, color = visc_colors[j], marker = viscmarker)
        end
        hlines!(axs[j], [f_bif ./ pa .+ f_pd], color = :gray10, linestyle = :dash, linewidth = 4)
    end
    vlines!(axs[j], log10.([dfdt_min, dfdt_max]), color = :gray20, linewidth = 3)

    axs[j].title = visc_num_labels[j]
    # j in (1,3,5) ? axs[j].xlabel = L"$\mathrm{log}_{10}$ slope (K)" : nothing
    j == 3 ? axs[j].xlabel = L"$\mathrm{log}_{10}$ rate of GMT anomaly, $\dot{f}$ $(\mathrm{K} \, \mathrm{yr}^{-1})$" : nothing
    axs[j].xticks = -4:1:0
    axs[j].xminorticks = -4:0.5:0
    axs[j].xminorticksvisible = true
    if region == "recovery"
        axs[j].yticks = 5.2:0.2:6.8 .+ f_pd
    elseif region == "wais"
        axs[j].yticks = 0.6:0.2:1.4 .+ f_pd
    end

    if region == "recovery"
        ylims!(axs[j], (9.3 ./ pa, 11.5 ./ pa) .+ f_pd)
        xlims!(axs[j], (log10(0.0003 ./ pa)-0.25, log10(0.3 ./ pa)+0.25))
    elseif region == "wais"
        ylims!(axs[j], (1.3 ./ pa, 2.7 ./ pa) .+ f_pd)
        xlims!(axs[j], (log10(3f-5 ./ pa)+0.25, log10(1 ./ pa)-0.15))
    end

    if 1 < j # < hr.n_visc_cases
        axs[j].yticksvisible = false
        axs[j].yticklabelsvisible = false
    end
end


for i in eachindex(f_hm_jet)
    scatter!(axs[k], log10.(hr.dfdt[k][idx_jet[i]] ./ pa),
        hr.f[k][idx_jet[i]] ./ pa .+ f_pd,
        markersize = ms2, color = jetmap[i])
    file1D = joinpath(hr.paths[k][idx_jet[i]], "yelmo1D.nc")
    lines!(ax_ramp1, ncread(file1D, "time") ./ 1e3,
        ncread(file1D, "hyst_f_now") ./ pa .+ f_pd,
        color = jetmap[i], linewidth = 2)
    lines!(ax_ramp2, ncread(file1D, "time") ./ 1e3,
        ncread(file1D, "hyst_f_now") ./ pa .+ f_pd,
        color = jetmap[i], linewidth = 2)
    lines!(ax_vol, ncread(file1D, "time")[1:100:end] ./ 1e3,
        ncread(file1D, "V_sle")[1:100:end],
        color = jetmap[i], linewidth = 2)
end

# The inset axis
X = ncread(joinpath(hr.paths[k][1], "yelmo2D.nc"), "x2D")
nx, ny = size(X)
XX, YY = ndgrid(1:nx, 1:ny)
dcrop = 10
ii = dcrop+1:nx-dcrop
jj = dcrop+1:ny-dcrop
valign_inset = region == "recovery" ? 0.15 : 0.85
inset_axs = [
    Axis(fig[2, i],
    # aspect = DataAspect(),
    width=Relative(0.3),
    height=Relative(0.3),
    halign=0.97,
    valign=valign_inset,
    backgroundcolor=:white) for i in [2, 3]]
[hidedecorations!(inset_ax) for inset_ax in inset_axs]

for i in eachindex(f_hm)
    
    hidedecorations!(ax_hm[i])
    file = joinpath(hr.paths[k][idx[i]], "yelmo2D.nc")
    time = ncread(file, "time")
    nt = length(time)
    H_ice_ref = ncread(file, "H_ice", start = [1, 1, 1], count = [-1, -1, 1])[:, :, 1]
    H_ice = ncread(file, "H_ice", start = [i1, i2, nt], count = [d1, d2, 1])[:, :, 1]
    z_bed = ncread(file, "z_bed", start = [i1, i2, nt], count = [d1, d2, 1])[:, :, 1]
    f_grnd_ref = ncread(file, "f_grnd", start = [i1, i2, 1], count = [d1, d2, 1])[:, :, 1]
    f_grnd = ncread(file, "f_grnd", start = [i1, i2, nt], count = [d1, d2, 1])[:, :, 1]
    f_grnd_ref_glob = ncread(file, "f_grnd", start = [1, 1, 1], count = [-1, -1, 1])[:, :, 1]
    heatmap!(ax_hm[i], z_bed; cmaps["z_bed"]...)
    heatmap!(ax_hm[i], H_ice + z_bed; cmaps["z_srf"]...)

    scatter!(axs[k], log10.(hr.dfdt[k][idx[i]] ./ pa), hr.f[k][idx[i]] ./ pa .+ f_pd,
        markersize = ms2, color = framecolors[i])
    file1D = joinpath(hr.paths[k][idx[i]], "yelmo1D.nc")
    contour!(ax_hm[i], f_grnd_ref, color = :black, linewidth = 2, levels = [0.5])
    contour!(ax_hm[i], f_grnd, color = :red, linewidth = 2, levels = [0.5])
    lines!(ax_ramp1, ncread(file1D, "time") ./ 1e3,
        ncread(file1D, "hyst_f_now") ./ pa .+ f_pd,
        color = framecolors[i], linewidth = lw)
    lines!(ax_ramp2, ncread(file1D, "time") ./ 1e3,
        ncread(file1D, "hyst_f_now") ./ pa .+ f_pd,
        color = framecolors[i], linewidth = lw)
    lines!(ax_vol, ncread(file1D, "time") ./ 1e3, ncread(file1D, "V_sle"),
        color = framecolors[i], linewidth = lw)
    heatmap!(inset_axs[i], H_ice_ref .> 1e-8, colorrange = (1e-8, 1),
        colormap = cgrad([:white, :gray70]), lowclip = :white, highclip = :gray70)
    contour!(inset_axs[i], f_grnd_ref_glob, color = :black, linewidth = 2, levels = [0.5])
    contour!(inset_axs[i], ((i1 .<= XX .<= i1+d1) .&& (i2 .<= YY .<= i2+d2)),
        color = :black, linewidth = 2)
    xlims!(inset_axs[i], 10, 371)
    ylims!(inset_axs[i], 10, 371)
    # text!(inset_axs[i], 10, 5, font = :bold, color = :black, fontsize = fs-4,
    #     text = "t=0 kyr")
        #text = "t = $(Int(round(time[nt] / 1e3, digits = 0))) kyr")
    ax_hm[i].leftspinecolor = framecolors[i]
    ax_hm[i].rightspinecolor = framecolors[i]
    ax_hm[i].topspinecolor = framecolors[i]
    ax_hm[i].bottomspinecolor = framecolors[i]
    ax_hm[i].spinewidth = 5
end

contour!(ax_hm[1], bbox[i1:i1+d1, i2:i2+d2], color = :darkorange, levels = [0.5], linewidth = 3)

if region == "recovery"
    text!(ax_ramp2, 35, 5 ./ pa .+ f_pd, text = "a", color = :black, font = :bold)
    text!(ax_vol, 37, 53, text = "b", color = :black, font = :bold)
    text!(ax_hm[1], 2, 102, text = "c", color = :white, font = :bold)
    text!(ax_hm[2], 2, 102, text = "d", color = :white, font = :bold)
    text!(ax_hm[1], 8, 103, text = "t=45 kyr", color = :white)
    text!(ax_hm[2], 8, 103, text = "t=45 kyr", color = :white)
elseif region == "wais"
    text!(ax_ramp1, 2.0, 1.1 ./ pa .+ f_pd, text = "a", color = :black, font = :bold)
    text!(ax_vol, 27, 57, text = "b", color = :black, font = :bold)
    text!(ax_hm[1], 2, 2, text = "c", color = :white, font = :bold)
    text!(ax_hm[2], 2, 2, text = "d", color = :white, font = :bold)
    text!(ax_hm[1], 8, 3, text = "t=30 kyr", color = :white)
    text!(ax_hm[2], 8, 3, text = "t=30 kyr", color = :white)
end

Colorbar(fig[1, 2], label = "Bed elevation (km)", vertical = false,
    width = Relative(rw_cbar), flipaxis = true, ticks = latexifyticks(-6:2:4, 1e3),
    halign = :left; cmaps["z_bed"]...)
Colorbar(fig[1, 2:3], label = "Ice surface elevation (km)", vertical = false,
    width = Relative(rw_cbar/2), flipaxis = true, ticks = (vcat(1, 1e3:1e3:4e3),
    string.(0:4)); cmaps["z_srf"]...)
Colorbar(fig[1, 3], label = L"$V_\mathrm{AIS}$ (m SLE)", vertical = false,
    width = Relative(rw_cbar), flipaxis = true, halign = :right; cmap...)

colgap!(ra, 0)
rowgap!(ra, 5)
rowsize!(ra, 1, 170)
rowsize!(ra, 2, 170)
colgap!(ha, 5)

colgap!(fig.layout, 10)
colsize!(fig.layout, 1, 200)
rowgap!(fig.layout, 1, -65)
rowgap!(fig.layout, 2, -30)
rowsize!(fig.layout, 2, 400)

# colgap!(fig.layout, 10)
# rowgap!(fig.layout, 10)

# rowgap!(fig.layout, 1, -60)
# rowgap!(fig.layout, 2, -40)
# rowgap!(fig.layout, 3, -30)
# rowsize!(fig.layout, 1, 10)

save(plotsdir("16km/rtip/ramp-heatmap-$region.png"), fig)
save(plotsdir("16km/rtip/ramp-heatmap-$region.pdf"), fig)


#####################################################################################
# Plot RSL for RSB and R-tip gap for other basins
#####################################################################################


function get_grline(f_grnd) 
    grline_cropped = (abs.(diff(f_grnd, dims = 1)) .> 0)[:, 2:end] .||
        (abs.(diff(f_grnd, dims = 2)) .> 0)[2:end, :]
    grline = zeros(Bool, size(f_grnd))
    grline[2:end, 2:end] .= grline_cropped
    return grline
end

function timeseries_from_var2D(t, X, mask::BitMatrix)
    x = zeros(eltype(X), length(t))
    for i in 1:length(t)
        x[i] = mean(X[mask, i])
    end
    return t, x
end

function timeseries_from_var2D(t, X, mask::Array{Bool, 3})
    x = zeros(eltype(X), length(t))
    for i in 1:length(t)
        x[i] = mean(X[mask[:, :, i], i])
    end
    return t, x
end

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

function Base.sort!(sr::StepRtip{T}) where T
    for i in 1:sr.n_visc_cases
        idx = sortperm(sr.f[i])
        sr.f[i] = sr.f[i][idx]
        sr.V[i] = sr.V[i][idx]
        sr.t_end[i] = sr.t_end[i][idx]
    end
    return sr
end

sr = StepRtip{Float32}(visc_cases)
prefix = "output/ais/ramps/16km"
dirs = [
    datadir("$prefix/steps-sigmarange"),
]
for dir in dirs
    aggregate_xp!(sr, dir)
end
sort!(sr)

basins = ["WAIS", "WSB low lat", "RSB", "WSB high lat", "ASB"]
# f_bif = [1.25, 4.3, 6.0, 7.11, 7.89]
# f_bif = [1.3, 4.4, 6.0, 7.11, 7.89]
f_bif_eq = [2.5, 5.6, 7.2, 8.3, 9.1]
f_bif_qeq = [2.45, 5.9, 7.3, 8.35, 9.15]
f_bif = f_bif_eq
V_bif = [55.5, 47.5, 34.5, 23.5, 13.5]
# V_bif = [55, 45, 34, 23, 13]

i_rtip = [[findfirst(V .< V_bif[i]) for i in eachindex(f_bif)] for V in sr.V]
f_rtip = [sr.f[i][i_rtip[i]] for i in eachindex(i_rtip)]
# f_bif = f_rtip[1]
@show i_rtip

plot_rsl_now = false
nrows = plot_rsl_now ? 3 : 2
ncols = 2
time_rsl = ncread(file_rsl_ref, "time")
nt_rsl = length(time_rsl)

slines_opts = (linewidth = lw, markersize = ms2)
set_theme!(theme_latexfonts())
fig2 = Figure(size=(900, 750), fontsize = fs+4)
axs = [Axis(fig2[i, j], aspect = AxisAspect(1.2)) for i in 1:nrows, j in 1:ncols]
axs[1, 1].title = "Case study: Recovery basin"
axs[1, 2].title = "R-tipping: all basins"

ssp_colors = [:darkblue, :darkorange, :red, :darkred]
xpos = [3.5, 3.5, 3.5, 1]
for i in eachindex(ssps)
    hlines!(axs[1, 2], maximum(ssps[i][:, 2]), color = ssp_colors[i], linewidth = 2)
    text!(axs[1, 2], xpos[i], maximum(ssps[i][:, 2]), text = ssp_labels[i],
        color = ssp_colors[i])
end

for i in eachindex(basins)
    scatterlines!(axs[1, 2], f_rtip[i] .+ f_pd, label = visc_labels[i],
        color = visc_colors[i]; slines_opts...)
end
scatterlines!(axs[1, 2], f_bif, label = "Bifurcation", color = :black,
    linestyle = :dash, marker = :cross; slines_opts...)
hlines!(axs[2, 2], 0, color = :gray70, linestyle = :dash, linewidth = lw)
for i in eachindex(basins)
    scatterlines!(axs[2, 2], f_rtip[i] .+ f_pd .- f_bif, label = visc_labels[i],
        color = visc_colors[i]; slines_opts...)
end

f_grnd_ref = ncslice(file_rsl_ref, "f_grnd", nt_rsl)
grline_ref = get_grline(f_grnd_ref) .& bbox

for i in eachindex(visc_cases)
    file_rsl_now = joinpath(hr.paths[i][rsl_index], "yelmo2D.nc")
    t = ncread(file_rsl_now, "time")
    Rsl = ncread(file_rsl_now, "z_sl") .- ncread(file_rsl_now, "z_bed")
    t, rsl_ref = timeseries_from_var2D(t, Rsl, grline_ref)

    if plot_rsl_now
        f_grnd_3D = ncread(file_rsl_now, "f_grnd")
        grline_3D = zeros(Bool, size(f_grnd_3D))
        for k in 1:nt_rsl
            grline_3D[:, :, k] .= get_grline(f_grnd_3D[:, :, k]) .& bbox
        end
        _, rsl_now = timeseries_from_var2D(t, Rsl, grline_3D)
        lines!(axs[3], t ./ 1f3, rsl_now, color = visc_colors[i],
            label = visc_labels[i], linewidth = lw)
    end

    t1D = ncread(joinpath(hr.paths[i][rsl_index], "yelmo1D.nc"), "time")
    V = ncread(joinpath(hr.paths[i][rsl_index], "yelmo1D.nc"), "V_sle")

    lines!(axs[1, 1], t ./ 1f3, rsl_ref, color = visc_colors[i],
        label = visc_labels[i], linewidth = lw)
    lines!(axs[2, 1], t1D ./ 1f3, V, color = visc_colors[i],
        label = visc_labels[i], linewidth = lw)
end

axs[1, 1].xticks = 0:5:25
axs[2, 1].xticks = 0:5:25
axs[1, 2].yaxisposition = :right
axs[1, 2].ylabel = "Critical forcing (K)"
axs[1, 2].ygridvisible = false
axs[1, 2].xgridvisible = false
axs[1, 2].xticklabelsvisible = false
axs[1, 2].yticks = 2:2:8

axs[2, 2].yaxisposition = :right
axs[2, 2].xticks = (eachindex(basins), basins)
axs[2, 2].xticklabelrotation = pi/6
axs[2, 2].xticklabelpad = 5
axs[2, 2].ylabel = "R-tipping gap (K)"
axs[2, 2].yticks = -1.6:0.2:0.1
axs[2, 2].xgridvisible = false
axs[2, 2].ygridvisible = false
ylims!(axs[1, 2], (1.8, 9.5))
ylims!(axs[2, 2], (-1.0, 0.2))
xlims!(axs[1, 2], (0.8, 5.2))
xlims!(axs[2, 2], (0.8, 5.2))

if plot_rsl_now
    axs[2].xticklabelsvisible = false
    axs[3].xticklabelsvisible = true
    axs[3].ylabel = "RSL tracked (m)"
    axs[3].xlabel = "Time (kyr)"
    xlims!(axs[3], (0, 30))
else
    axs[2, 1].xticklabelsvisible = true
    axs[2, 1].xlabel = "Time (kyr)"
end

axs[1, 1].xticklabelsvisible = false
axs[1, 1].ylabel = "Mean relative SL (m)"
axs[2, 1].ylabel = "AIS volume (m SLE)"
xlims!(axs[1, 1], (0, 30))
xlims!(axs[2, 1], (0, 30))
Legend(fig2[0, 1:2], axs[1, 2], nbanks = 6)
rowsize!(fig2.layout, 0, 20)

text!(axs[1, 1], 27, 300, text="a", font=:bold)
text!(axs[2, 1], 27, 55, text="b", font=:bold)
text!(axs[1, 2], 0.95, 8.7, text="c", font=:bold)
text!(axs[2, 2], 0.95, 0.05, text="d", font=:bold)

axs[1, 1].xticksvisible = false
axs[1, 2].xticksvisible = false
rowgap!(fig2.layout, 5)
rowgap!(fig2.layout, 1, 10)
colgap!(fig2.layout, 5)
save(plotsdir("16km/rtip/rsl-$region.png"), fig2)
save(plotsdir("16km/rtip/rsl-$region.pdf"), fig2)