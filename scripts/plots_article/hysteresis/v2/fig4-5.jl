include("../../../intro.jl")

T = Float32
regrowth_dir = datadir("output/ais/v2/hyster/regrowth")
xps = [
    "$regrowth_dir/aqef/dpr",
    "$regrowth_dir/aqef/upl",
    "$regrowth_dir/aqef/refnomslow",
]
aqef = AQEFResults(T, xps)
# eql = EquilResults(T, "$regrowth_dir/equil")

xp_labels = [
    "DPR",
    "UPL",
    "REF",
]
lw1, lw2 = 3, 6
lws = [lw1, lw1, lw2]
cycling_colors = [
    xpcolors["DPR"],
    xpcolors["UPL"],
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
cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 18
fig4 = Figure(size=(800, 600), fontsize = 24)
ax = Axis(fig4[1, 1])
ax.xticks = -2:1:12
ax.xminorticks = -2:0.2:12
ax.yticks = 0:10:60
ax.yminorticks = IntervalsBetween(10)
ax.xminorgridvisible = true
ax.yminorgridvisible = true
ax.xlabel = L"GMT anomaly $f$ (K)"
ax.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(ax, 0, 60)
xlims!(ax, 1, 11)

bifs = [
    (9.0, "East Gamburtsev", 30, 15.5),
    (8.657, "Gamburtsev-Droning Maud", 8, 24),
    (8.462, "Gamburtsev-Transantarctic", 33, 24),
    (8.153, "Droning Maud", 12, 13),
    (6.457, "Recovery-Gamburtsev-Aurora", 15, 25.7),
    (4.853, "Transantarctic", 20, 13),
    (3.8, "Pensacola-Pole, Amery & Recovery", 27, 30),
    # (3.738, "Recovery & Amery", 50, 30),
    (3.21, "Recovery", 12, 8.5),
    (2.836, "Pensacola-Pole", 40, 13.2),
    (2.658, "Aurora", 17, 6.8),
    (2.573, "Aurora-Wilkes", 45, 13),
    (2.424, "Wilkes", 25, 6.5),
    (1.174, "Wilkes", 30, 6.5),
]

f_bif = [bifs[i][1] for i in eachindex(bifs)]
vlines!(ax, f_bif, color = :gray60, alpha = 0.5, linewidth = 5)

for i in eachindex(bifs)
    rectangle!(ax, (bifs[i][1] - 0.13, bifs[i][3] - 0.4), 0.3, bifs[i][4])
    text!(ax, bifs[i][1] + 0.15, bifs[i][3], text = bifs[i][2], rotation = π/2, fontsize = 16)
end
vlines!(ax, 1f6, alpha = 0.9, color = :gray60, linewidth = 6, label = "Bifurcation")

s = 200
alpha = 1
for k in 1:aqef.n_xps
    lines!(ax, aqef.f[k][1:s:end] ./ polar_amplification .+ f2015,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]), alpha = alpha)
end
# scatter!(axs[1, 1], eql.f ./ polar_amplification .+ f2015, eql.V_sle;
#     color = :black, label = "EQL", markersize = ms1)
axislegend(ax, position = :rt)
fig4
save(plotsdir("v2/hysteresis/fig4.png"), fig4)
save(plotsdir("v2/hysteresis/fig4.pdf"), fig4)


##############################################
# Fig 5
###############################################

fig5 = Figure(size=(1600, 1200), fontsize = 24)
nrows, ncols = 3, 4
axs = [Axis(fig5[i+1, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]
f_map = f_bif[vcat(2, 4, 7, 8, 10, 12)]
f_map = vcat(f_map .+ 0.05, f_map .- 0.05)
sort!(f_map, rev = true)
forcing_frames = permutedims(reshape(f_map, ncols, nrows))
state_labels = ["a" "b" "c" "d";
    "e" "f" "g" "h";
    "i" "j" "k" "l"]

file2D = joinpath(aqef.xps[end], "0", "yelmo2D.nc")
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))

ii = cropx+1:nx-cropx
jj = cropy+1:ny-cropy
XX = X[ii, jj]
YY = Y[ii, jj]

# xlims_frames = permutedims(reshape([nothing, (500, 1500), (500, 1500),
#     nothing, nothing, nothing, nothing, nothing, nothing,
#     nothing, nothing, nothing], ncols, nrows))

# ylims_frames = permutedims(reshape([nothing, (300, 1300), (300, 1300),
#     nothing, nothing, nothing, nothing, nothing, nothing,
#     nothing, nothing, nothing], ncols, nrows))

var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    forcing = forcing_frames[i, j]
    file2D_plot = file2D
    plot_index = xp_idx

    i3 = argmin( ( aqef.f[plot_index] ./ polar_amplification .+ f2015 .- forcing) .^ 2)
    f_eq = aqef.f[plot_index][i3] ./ polar_amplification .+ f2015
    V_eq = aqef.V_sle[plot_index][i3]
    frame_index = argmin(abs.(aqef.t_2D[plot_index] .- aqef.t_1D[plot_index][i3]))
    z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D_plot,
        var_names_2D, frame_index)
   
    @show i3, f_eq, V_eq

    hidedecorations!(axs[i, j])
    heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
    heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
    contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
        color = :red, linewidth = 2)
    # if xlims_frames[i, j] !== nothing
    #     contour!(axs[i, j], xc, yc, (xlims_frames[i, j][1] .< X .< xlims_frames[i, j][2]) .&
    #         (ylims_frames[i, j][1] .< Y .< ylims_frames[i, j][2]), levels = [0.5],
    #         color = :darkred, linewidth = 3)
    # end
    text!(axs[i, j], -2500, -2450, color = :black, font = :bold,
        text="("*state_labels[i, j]*") f = $(round(forcing, digits=2)) K", fontsize = 24)
    xlims!(axs[i, j], extrema(XX))
    ylims!(axs[i, j], extrema(YY))
end


relwidth = 0.4
Colorbar(fig5[1, 1:2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig5[1, 3:4], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :darkred, linewidth = 2)
Legend(fig5[1, 2:3], [elem_1, elem_2], ["Grounding line", "Highlighted region"], valign = :bottom)

rowsize_base = 350
rowgap!(fig5.layout, 5)
colgap!(fig5.layout, 5)
rowgap!(fig5.layout, 1, 30)
rowsize!(fig5.layout, 1, 1)
rowsize!(fig5.layout, 2, rowsize_base)
rowsize!(fig5.layout, 3, rowsize_base)
rowsize!(fig5.layout, 4, rowsize_base)
colsize!(fig5.layout, 1, rowsize_base*aratio)
colsize!(fig5.layout, 2, rowsize_base*aratio)
colsize!(fig5.layout, 3, rowsize_base*aratio)
colsize!(fig5.layout, 4, rowsize_base*aratio)
fig5
save(plotsdir("v2/hysteresis/fig5.png"), fig5)
save(plotsdir("v2/hysteresis/fig5.pdf"), fig5)