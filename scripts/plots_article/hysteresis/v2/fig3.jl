include("../../intro.jl")

T = Float32
visc_type = "normvisc"    # "equil" or "aqef"
regrowth_dir = datadir("output/ais/hyster/16km/regrowth")
xps = [
    "$regrowth_dir/aqef/pmpt-$visc_type-dpr",
    "$regrowth_dir/aqef/pmpt-$visc_type-upl",
    "$regrowth_dir/aqef/pmpt-$visc_type-fastnormforcing-restarted",
    "$regrowth_dir/aqef/pmpt-$visc_type-fastnormforcing",
]
aqef = AQEFResults(T, xps)
eqldir = datadir("output/ais/hyster/16km/regrowth/equil")
eql = EquilResults(T, eqldir)

xp_labels = [
    "DPR",
    "UPL",
    nothing,
    "REF",
]
lw1, lw2 = 3, 6
lws = [lw1, lw1, lw2, lw2]
cycling_colors = [
    xpcolors["DPR"],
    xpcolors["UPL"],
    xpcolors["REF"],
    xpcolors["REF"],
]

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
f2015 = 1.2
heatmap_frames = "aqef"    # "equil" or "aqef"
xp_idx = aqef.n_xps
f_ref = aqef.f[end] ./ polar_amplification
stiching_forcing = 2
stiching_idx = findfirst(f_ref .<= stiching_forcing)
cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 18
nrows, ncols = 3, 4
forcing_frames = permutedims(reshape([nothing, 7.6, 7.2, 6.8, 3.8, 3, 2.2, 1.8, 1.6, 1, 0.2, -2],
    ncols, nrows))
state_labels = ["a" "b" "c" "d";
    "e" "f" "g" "h";
    "i" "j" "k" "l"]
shade = [(1.3, 1.4), (2.6, 2.7), (2.9, 3.0), (3.25, 3.35),
    (3.6, 3.7), (4.1, 4.2), (4.8, 4.9), (8.05, 8.15), (8.65, 8.75)]

fig = Figure(size=(1400, 1050), fontsize = 24)
axs = [Axis(fig[i+1, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]
shade_transitions!([axs[1, 1]], shade, 0.0, 0.2, :gray70)
vlines!(axs[1, 1], 1f6, color = dgray, linewidth = 5, alpha = 0.6, label = "Bifurcation")

axs[1, 1].xticks = -2:2:12
axs[1, 1].xminorticks = -2:0.2:12
axs[1, 1].yticks = 0:10:60
axs[1, 1].yminorticks = IntervalsBetween(10)
axs[1, 1].xminorgridvisible = true
axs[1, 1].yminorgridvisible = true
axs[1, 1].xaxisposition = :top
axs[1, 1].xlabel = L"GMT anomaly $f$ (K)"
axs[1, 1].ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(axs[1, 1], 0, 60)
xlims!(axs[1, 1], -1, 11)
fig

file2D_restarted = joinpath(aqef.xps[xp_idx - 1], "0", "yelmo2D.nc")
file2D = joinpath(aqef.xps[xp_idx], "0", "yelmo2D.nc")
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))

ii = cropx+1:nx-cropx
jj = cropy+1:ny-cropy
XX = X[ii, jj]
YY = Y[ii, jj]

text_offsets = permutedims(reshape([
    nothing,         # 1
    (-25, -5),         # 2
    (0, -35),      # 3
    (-27, -3),     # 4
    (-10, -35),       # 5
    (-9, -37),       # 6
    (-8, -35),     # 7
    (10, -15),       # 8
    (10, -15),       # 9
    (-27, 0),     # 10
    (-5, -37),       # 11
    (-10, 0),       # 12
], ncols, nrows))

xlims_frames = permutedims(reshape([nothing, (500, 1500), (500, 1500),
    nothing, nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing], ncols, nrows))

ylims_frames = permutedims(reshape([nothing, (300, 1300), (300, 1300),
    nothing, nothing, nothing, nothing, nothing, nothing,
    nothing, nothing, nothing], ncols, nrows))

var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    if i > 1 || j > 1
        forcing = forcing_frames[i, j]
        if heatmap_frames == "aqef" || (i == 1 && j == 2)
            
            if forcing > stiching_forcing
                file2D_plot = file2D
                plot_index = xp_idx
            else
                file2D_plot = file2D_restarted
                plot_index = xp_idx - 1
            end

            i3 = argmin( ( aqef.f[plot_index] ./ polar_amplification .-
                forcing) .^ 2)
            f_eq = aqef.f[plot_index][i3] ./ polar_amplification
            V_eq = aqef.V_sle[plot_index][i3]
            frame_index = argmin(abs.(aqef.t_2D[plot_index] .- aqef.t_1D[plot_index][i3]))
            z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D_plot,
                var_names_2D, frame_index)
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

        scatter!(axs[1, 1], f_eq .+ f2015, V_eq, color = :red, markersize = ms2)
        text!(axs[1, 1], f_eq .+ f2015, V_eq, text = state_labels[i, j],
            color = :red, font = :bold, fontsize = 30, offset = text_offsets[i, j])

        hidedecorations!(axs[i, j])
        heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
        heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
        contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
            color = :red, linewidth = 2)
        if xlims_frames[i, j] !== nothing
            contour!(axs[i, j], xc, yc, (xlims_frames[i, j][1] .< X .< xlims_frames[i, j][2]) .&
                (ylims_frames[i, j][1] .< Y .< ylims_frames[i, j][2]), levels = [0.5],
                color = :darkred, linewidth = 3)
        end
        text!(axs[i, j], -2500, -2450, color = :white, font = :bold,
            text="("*state_labels[i, j]*")", fontsize = 30)
        xlims!(axs[i, j], extrema(XX))
        ylims!(axs[i, j], extrema(YY))
    end
end
s = 200
alpha = 1
for k in 1:aqef.n_xps
    if k < aqef.n_xps
        lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification .+ f2015,
            aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
            color = lcolor(cycling_colors[k]), alpha = alpha)
    else
        lines!(axs[1, 1], aqef.f[k][1:s:stiching_idx] ./ polar_amplification .+ f2015,
            aqef.V_sle[k][1:s:stiching_idx], linewidth = lws[k], label = xp_labels[k],
            color = lcolor(cycling_colors[k]), alpha = alpha)
    end
end

# f2015 = 1.12
scatter!(axs[1, 1], eql.f ./ polar_amplification .+ f2015, eql.V_sle;
    color = :black, label = "EQL", markersize = ms1)

text!(axs[1, 1], 1, 1, font = :bold, text = "(a)", color = :black, fontsize = 30)
axislegend(axs[1, 1], position = :lt, nbanks = 1)
relwidth = 0.8
Colorbar(fig[1, 2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[1, 3], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :darkred, linewidth = 2)
Legend(fig[1, 4], [elem_1, elem_2], ["Grounding line", "Highlighted region"])

axs[1, 1].xreversed = true
rowsize_base = 300
rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -20)
rowsize!(fig.layout, 1, 1)
rowsize!(fig.layout, 2, rowsize_base)
rowsize!(fig.layout, 3, rowsize_base)
rowsize!(fig.layout, 4, rowsize_base)
colsize!(fig.layout, 1, rowsize_base*aratio)
colsize!(fig.layout, 2, rowsize_base*aratio)
colsize!(fig.layout, 3, rowsize_base*aratio)
colsize!(fig.layout, 4, rowsize_base*aratio)
fig
save(plotsdir("16km/hysteresis/fig3.png"), fig)
save(plotsdir("16km/hysteresis/fig3.pdf"), fig)