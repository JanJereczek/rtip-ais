include("../../../intro.jl")

T = Float32
regrowth_dir = datadir("output/ais/v2/hyster/regrowth")
xps = [
    "$regrowth_dir/aqef/dpr",
    "$regrowth_dir/aqef/upl",
    "$regrowth_dir/aqef/minvisc/refnomslow-restarted",
    "$regrowth_dir/aqef/refnomslow",
]
aqef = AQEFResults(T, xps)
eql = EquilResults(T, "$regrowth_dir/equil/refnomslow")

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
f2020 = 1.2
heatmap_frames = "aqef"    # "equil" or "aqef"
xp_idx = aqef.n_xps
f_ref = aqef.f[end] ./ polar_amplification
cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 18
fig4 = Figure(size=(800, 800), fontsize = 24)
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
xlims!(ax, 0, 10)

bifs = [
    (9.0, "East Gamburtsev", 5, 19),
    (8.657, "Gamburtsev-Queen Maud", 8, 27),
    # (8.462, "Gamburtsev-Transantarctic", 33, 24),
    (8.153, "Queen Maud", 12, 13.8),
    (6.457, "Recovery-Gamburtsev-Aurora", 16, 31),
    (4.853, "Transantarctic", 20, 15.5),
    (3.8, "Pensacola-Pole, Amery & Recovery", 23, 36.5),
    (3.5, "Vostok", 10, 8),
    (3.21, "Recovery", 37, 10),
    (2.836, "Pensacola-Pole", 40, 16),
    (2.658, "Aurora", 17, 7.8),
    (2.573, "Aurora-Wilkes", 45, 15.5),
    (2.424, "Wilkes", 25, 7.5),
    (1.174, "Wilkes", 35, 7.5),
    (0.9, "Wilkes & Siple Coast", 5, 22),
    (0.45, "Wilkes", 40, 7.5),
]

f_bif = [bifs[i][1] for i in eachindex(bifs)]
vlines!(ax, f_bif, color = :gray60, alpha = 0.5, linewidth = 5)

for i in eachindex(bifs)
    rectangle!(ax, (bifs[i][1] - 0.13, bifs[i][3] - 0.4), 0.3, 0.7*bifs[i][4])
    text!(ax, bifs[i][1] - 0.12, bifs[i][3], text = bifs[i][2], rotation = π/2, fontsize = 18)
end
vlines!(ax, 1f6, alpha = 0.9, color = :gray60, linewidth = 6, label = "Bifurcation")

s = 200
alpha = 1
for k in 1:aqef.n_xps
    lines!(ax, aqef.f[k][1:s:end] ./ polar_amplification .+ f2020,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]), alpha = alpha)
end

idx = sortperm(eql.f)
f = eql.f[idx]
V_sle = eql.V_sle[idx]
V_sle = filter_outliers(V_sle, mode = :low)
scatter!(ax, f ./ polar_amplification .+ f2020, V_sle;
    color = :cornflowerblue, label = "EQL", markersize = ms1)
axislegend(ax, position = :lt)
# Legend(fig4[0, 1], ax, nbanks = 5, framevisible = false)

# Compute mean diff between eql and ref
idx = sortperm(aqef.f[end])
aqef_itp = linear_interpolation(aqef.f[end][idx] ./ polar_amplification .+ f2020, aqef.V_sle[end][idx])
aqef_eql = aqef_itp.(f ./ polar_amplification .+ f2020)
mean_diff = mean(abs.(aqef_eql[f .< 16] .- V_sle[f .< 16])) # exclude ranges where volue = 0 to avoid underestimating error
println("Mean difference REF vs EQL: $(round(mean_diff, digits=2)) m SLE")
lines(f[f .< 16], aqef_eql[f .< 16] .- V_sle[f .< 16])
ax.xreversed = true

# rowsize!(fig4.layout, 0, 20)
colsize!(fig4.layout, 1, 700)
# axislegend(ax, position = :ct)
fig4
save(plotsdir("v2/hysteresis/fig4.png"), fig4)
save(plotsdir("v2/hysteresis/fig4.pdf"), fig4)


##############################################
# Fig 5
###############################################

fig5 = Figure(size=(1600, 1200), fontsize = 24)
nrows, ncols = 3, 4
axs = [Axis(fig5[i+1, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]
f_map = f_bif[vcat(2, 3, 6, 8, 10, 12)]
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

xl = [
    [(400, 1500), (800, 1800)],    # Gamburtsev-Queen Maud
    [(-300, 500), (1400, 2200)],      # Queen Maud
    [(-600, 400), (1200, 2000), (-400, 400)],    # Pensacola-Pole
    [(-600, 500)],    # Recovery
    [(1200, 2500)],   # Aurora
    [(500, 1500)],    # Wilkes
]
yl = [
    [(300, 1400), (-400, 400)],    # Gamburtsev-Queen Maud
    [(1200, 2100), (800, 1800)],   # Queen Maud
    [(-100, 500), (300, 1100), (1000, 1600)],       # Pensacola-Pole
    [(400, 1500)],    # Recovery
    [(-1200, -200)],  # Aurora
    [(-2000, -800)],  # Wilkes
]
x_highlight = permutedims(reshape(repeat(xl, inner = 2), ncols, nrows), (2, 1))
y_highlight = permutedims(reshape(repeat(yl, inner = 2), ncols, nrows), (2, 1))
var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    forcing = forcing_frames[i, j]
    file2D_plot = file2D
    plot_index = xp_idx

    i3 = argmin( ( aqef.f[plot_index] ./ polar_amplification .+ f2020 .- forcing) .^ 2)
    f_eq = aqef.f[plot_index][i3] ./ polar_amplification .+ f2020
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
    for k in eachindex(x_highlight[i, j])
        xl = x_highlight[i, j][k]
        yl = y_highlight[i, j][k]
        contour!(axs[i, j], xc, yc, (xl[1] .< X .< xl[2]) .&
            (yl[1] .< Y .< yl[2]), levels = [0.5],
            color = :darkred, linewidth = 3)
    end
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
Legend(fig5[1, 2:3], [elem_1, elem_2], ["Grounding line", "Highlighted regions"], valign = :bottom)

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