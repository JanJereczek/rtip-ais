include("../intro.jl")


T = Float32
visc_type = "normvisc"    # "equil" or "aqef"
xps = [
    "pmpt-$visc_type-ocnforcing",
    "pmpt-$visc_type-atmforcing",
    "pmpt-$visc_type-normforcing",
]
xp_labels = [
    "OCN",
    "ATM",
    "REF",
]
lws = [3, 3, 5]
aqef = AQEFResults(T, datadir("output/ais/hyster/16km/retreat/aqef"), xps)
eqldir = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
eql = EquilResults(T, eqldir)

############################################################################
# A simple mb plot for the appendix
############################################################################
plot_mb = false

if plot_mb
    lw = 3
    set_theme!(theme_latexfonts())
    fig = Figure(size = (800, 500), fontsize = 24)

    ax2 = Axis(fig[1, 1], yaxisposition = :right)
    hidexdecorations!(ax2)
    dvdt = diff(aqef.V_sle[1]) ./ diff(aqef.t_1D[1])
    s = 10
    lines!(ax2, aqef.t_1D[1][2:s:end] ./ 1f3, dvdt[1:s:end] .* 1f3, label = "Volume change",
        linewidth = lw, color = :gray70)
    axislegend(ax2, position = :rb, labelsize = 20)
    xlims!(ax2, 0, 200)
    ylims!(ax2, -2, 2)
    ax2.ylabel = "Volume change (mmSLE/yr)"
    ax2.yticks = -2:1:2

    ax = Axis(fig[1, 1])
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.bmb[1], label = "basal", linewidth = lw)
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.smb[1], label = "surface", linewidth = lw)
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.cmb[1], label = "calving", linewidth = lw)
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.bmb[1] .+ aqef.smb[1] .+ aqef.cmb[1], 
    linewidth = lw, label = "total")
    axislegend(ax, position = :rt, labelsize = 20)
    xlims!(ax, 0, 200)
    ylims!(ax, -0.1, 0.1)
    ax.xlabel = "Time (kyr)"
    ax.ylabel = "Mass balance (m/yr)"
    ax.yticks = -0.1:0.05:0.1

    save(plotsdir("16km/hysteresis/fluxes.png"), fig)
end

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
heatmap_frames = "equil"    # "equil" or "aqef"
xp_idx = aqef.n_xps

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 15
nrows, ncols = 3, 4
forcing_frames = reshape([0, 0, 1.15, 1.35, 4.2, 4.4, 5.8, 6.0, 6.7, 7.1, 7.6, 7.8], ncols, nrows)'
state_labels = latexify.(reshape(0:11, ncols, nrows)')
fig = Figure(size=(1700, 1320), fontsize = 24)
axs = [Axis(fig[i+1, j], aspect = AxisAspect(1)) for i in 1:nrows, j in 1:ncols]
s = 50
for k in 1:aqef.n_xps
    lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k])
end
scatterlines!(axs[1, 1], eql.f ./ polar_amplification, eql.V_sle;
    linewidth = lws[xp_idx], color = :black, label = "EQL", markersize = ms1)
axs[1, 1].xticks = 0:2:12
axs[1, 1].xminorticks = 0:0.2:12
axs[1, 1].yticks = 0:10:60
axs[1, 1].yminorticks = IntervalsBetween(10)
axs[1, 1].xminorgridvisible = true
axs[1, 1].yminorgridvisible = true
axs[1, 1].xaxisposition = :top
axs[1, 1].xlabel = L"GMT anomaly $f$ (K)"
axs[1, 1].ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(axs[1, 1], 0, 60)
xlims!(axs[1, 1], 0, 11)
fig

plot_ssp = false
if plot_ssp

    hist = readdlm(datadir("processed/SSP/History.csv"), ',')
    f2014 = hist[end, 2]
    ssp1 = readdlm(datadir("processed/SSP/SSP1.csv"), ',') .- f2014
    ssp2 = readdlm(datadir("processed/SSP/SSP2.csv"), ',') .- f2014
    ssp3 = readdlm(datadir("processed/SSP/SSP3.csv"), ',') .- f2014
    # ssp4 = readdlm(datadir("processed/SSP/SSP4.csv"), ',')
    ssp5 = readdlm(datadir("processed/SSP/SSP5.csv"), ',') .- f2014
    ssp1_2100 = ssp1[end, 2]
    ssp2_2100 = ssp2[end, 2]
    ssp3_2100 = ssp3[end, 2]
    ssp5_2100 = ssp5[end, 2]

    line_opts = (linewidth = 3, linestyle = :dash)
    vlines!(axs[1, 1], [ssp1_2100], color = :darkblue, label = "SSP1-2100"; line_opts...)
    vlines!(axs[1, 1], [ssp2_2100], color = :lightblue, label = "SSP2-2100"; line_opts...)
    vlines!(axs[1, 1], [ssp3_2100], color = :orange, label = "SSP3-2100"; line_opts...)
    vlines!(axs[1, 1], [ssp5_2100], color = :darkred, label = "SSP5-2100"; line_opts...)
end

# equil_opts = (label = "Equilibrium", color = :black, linewidth = 4)
# scatterlines!(axs[1, 1], f_equil, V_equil; equil_opts...)

file2D = joinpath(aqef.dir, aqef.xps[xp_idx], "0", "yelmo2D.nc")
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))

crop = 20
ii = crop+1:nx-crop
jj = crop+1:ny-crop
XX = X[ii, jj]
YY = Y[ii, jj]

struct Transect
    i1::Int
    i2::Int
    j1::Int
    j2::Int
end
makeline(ax, tr::Transect; kwargs...) = lines!(ax, [tr.i1, tr.i2], [tr.j1, tr.j2]; kwargs...)
tr1 = Transect(69, 100, 140, 160)
tr2 = Transect(242, 240, 35, 70)
tr3 = Transect(148, 160, 247, 255)
tr4 = Transect(237, 220, 72, 120)
tr5 = Transect(320, 270, 128, 128)
transects_vec = [nothing, tr1, tr2, tr3, tr4, tr5]

struct PrecomputedTransect{T<:AbstractFloat}
    tr::Transect
    x::Vector{T}
    y::Vector{T}
    dist::Vector{T}
end

function line_on(tr::Nothing, X, Y)
    return nothing
end

function line_on(tr::Transect, X, Y)

    ii = [tr.i1, tr.i2]
    jj = [tr.j1, tr.j2]
    i1 = minimum(ii)
    i2 = maximum(ii)
    j1 = minimum(jj)
    j2 = maximum(jj)

    x1, x2 = X[tr.i1, tr.j1], X[tr.i2, tr.j2]
    y1, y2 = Y[tr.i1, tr.j1], Y[tr.i2, tr.j2]

    if abs(tr.i1 - tr.i2) > abs(tr.j1 - tr.j2)
        x = X[i1:i2, 1]
        m = (y2 - y1) / (x2 - x1)
        p = y1 - m * x1
        y = m .* x .+ p
    else
        y = Y[1, j1:j2]
        m = (x2 - x1) / (y2 - y1)
        p = x1 - m * y1
        x = m .* y .+ p
    end

    dist = sqrt.((x .- X[tr.i1, tr.j1]).^2 .+ (y .- Y[tr.i1, tr.j1]).^2)
    idx = sortperm(dist)
    
    return PrecomputedTransect(tr, x[idx], y[idx], dist[idx])
end

ntransects = length(transects_vec)
ptrs = [line_on(tr, XX, YY) for tr in transects_vec]
transects = permutedims(reshape(repeat(ptrs, inner = 2), ncols, nrows))

# var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
# velmap = (colormap = :inferno, colorrange = (1e2, 1e4), lowclip = :transparent)

text_offsets = permutedims(reshape([
    (0, 0),         # 0
    (7, -33),         # 1
    (7, -12),      # 2
    (-20, -32),     # 3
    (7, -10),       # 4
    (-17, -10),       # 5
    (7, -11),     # 6
    (-18, -20),       # 7
    (7, -15),       # 8
    (-22, -15),     # 9
    (5, -10),       # 10
    (10, -22),       # 11
], ncols, nrows))

var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    if i > 1 || j > 1
        forcing = forcing_frames[i, j]
        if heatmap_frames == "aqef" || (i == 1 && j == 2)
            i3 = findfirst(aqef.f[xp_idx] ./ polar_amplification .>= forcing)
            f_eq, V_eq = aqef.f[xp_idx][i3] ./ polar_amplification, aqef.V_sle[xp_idx][i3]
            frame_index = argmin(abs.(aqef.t_2D[xp_idx] .- aqef.t_1D[xp_idx][i3]))
            z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D,
                frame_index)
        else
            i3 = argmin((eql.f ./ polar_amplification .- forcing) .^ 2)
            file_2D_equil = datadir("$eqldir/$(string(i3))/0/yelmo2D.nc")
            @show file_2D_equil
            f_eq, V_eq = eql.f[i3] ./polar_amplification, eql.V_sle[i3]
            frame_index = length(ncread(file_2D_equil, "time"))
            z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file_2D_equil,
                var_names_2D, frame_index)
        end
        @show i3, f_eq, V_eq

        if mod(j, 2) == 1
            if i == 1 && j == 3
                vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.2, color = :gray, label = "Bifurcations")
            else
                vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.2, color = :gray)
            end
        end

        scatter!(axs[1, 1], f_eq, V_eq, color = :red, markersize = ms2)
        text!(axs[1, 1], f_eq, V_eq, text = state_labels[i, j],
            color = :grey10, fontsize = 30, font = :bold, offset = text_offsets[i, j])

        hidedecorations!(axs[i, j])
        heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
        heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
        contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
            color = :red, linewidth = 2)
        if transects[i, j] !== nothing
            lines!(axs[i, j], transects[i, j].x, transects[i, j].y; color = :orange, linewidth = 3)
        end
        text!(axs[i, j], -3000, -3000, color = :white, font = :bold,
            text=state_labels[i, j], fontsize = 30)
        xlims!(axs[i, j], extrema(XX))
        ylims!(axs[i, j], extrema(YY))
    end
end

axislegend(axs[1, 1], position = :lb, nbanks = 1)
relwidth = 0.8
Colorbar(fig[1, 2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[1, 3], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :orange, linewidth = 3)
Legend(fig[1, 4], [elem_1, elem_2], ["Grounding line", "Transects"])

rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -20)
rowsize!(fig.layout, 1, 1)
rowsize!(fig.layout, 2, 400)
rowsize!(fig.layout, 3, 400)
rowsize!(fig.layout, 4, 400)
colsize!(fig.layout, 1, 400)
colsize!(fig.layout, 2, 400)
colsize!(fig.layout, 3, 400)
colsize!(fig.layout, 4, 400)
save(plotsdir("16km/hysteresis/retreat.png"), fig)

transect_fig = Figure(size = (800, 2000), fontsize = 24)
axs = [Axis(transect_fig[i, 1], aspect = AxisAspect(2)) for i in 1:5]

bif_time_1 = [28, 90, 173, 243, 310] .* 1f3
bif_time_2 = bif_time_1 .+ 10f3

for i in 1:5
    @show i
    k1 = argmin(abs.(aqef.t_2D[xp_idx] .- bif_time_1[i]))
    k2 = argmin(abs.(aqef.t_2D[xp_idx] .- bif_time_2[i]))
    z_bed, z_srf, f_ice, f_grnd = load_netcdf_2D(file2D, ["z_bed", "z_srf", "f_ice", "f_grnd"], k1)
    
    z_bed_itp = linear_interpolation((xc, yc), z_bed)
    z_srf_itp = linear_interpolation((xc, yc), z_srf)
    f_ice_itp = linear_interpolation((xc, yc), f_ice)
    f_grnd_itp = linear_interpolation((xc, yc), f_grnd)

    ptr = ptrs[i+1]
    
    z_bed_vec = [z_bed_itp(ptr.x[j], ptr.y[j]) for j in 1:length(ptr.x)]
    z_srf_vec = [z_srf_itp(ptr.x[j], ptr.y[j]) for j in 1:length(ptr.x)]
    f_ice_vec = [f_ice_itp(ptr.x[j], ptr.y[j]) for j in 1:length(ptr.x)]
    f_grnd_vec = [f_grnd_itp(ptr.x[j], ptr.y[j]) for j in 1:length(ptr.x)]

    lines!(axs[i], ptr.dist, z_bed_vec, color = :black, linewidth = 3)
    lines!(axs[i], ptr.dist, z_srf_vec, color = :red, linewidth = 3)
    ylims!(axs[i], -1000, 1000)
end
transect_fig