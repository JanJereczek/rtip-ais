include("../../../intro.jl")

T = Float32
heatmap_frames = "aqef"    # "equil" or "aqef"
xps = [
    datadir("output/ais/v2/hyster/retreat/aqef/minvisc/refnomslow"),
]
xp_labels = [
    "REF",
]
aqef = AQEFResults(T, xps)
eqldir = datadir("output/ais/v2/hyster/retreat/equil/refnomslow")
eql = EquilResults(T, eqldir)
lws = [10]
cycling_colors = [xpcolors[l] for l in xp_labels]

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
xp_idx = aqef.n_xps
f2020 = 1.2

ms1, ms2 = 10, 15
s = 400

fn_ref = "/p/projects/megarun/ice_data/Antarctica/ANT-16KM/ANT-16KM_TOPO-RTOPO-2.0.1.nc"
mask_ref = ncread(fn_ref, "mask")
mask_ref[250:300, 150:200] .= 2

file2D = joinpath(aqef.xps[xp_idx], "0", "yelmo2D.nc")
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))
var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D, 1)

###############################
# Fig 1
###############################
alpha = 1
set_theme!(theme_latexfonts())
fig1 = Figure(size=(1600, 800), fontsize = 24)
ax = Axis(fig1[1, 1])
large_misi = [1.95, 5.87, 6.95, 8.65]
small_misi = [4.2, 6.21, 7.2, 7.6, 7.9, 9.2]
perimeter = [9.7, 10.05]
shadings = [large_misi, small_misi, perimeter]
colors = [:gray60, :gray60, :gray60]
labels = ["Large MISI", "Small MISI", "Perimeter"]
for i in eachindex(shadings)
    for j in eachindex(shadings[i])
        vlines!(ax, shadings[i][j], linewidth = 5, color = :gray60, alpha = 0.5)
    end
end
vlines!(ax, 1f6, alpha = 0.9, color = :gray60, linewidth = 6, label = "Bifurcation")

bifurcation_data = [
    ("WAIS", 6.5, 38),
    ("Amery", 7.0, 38),
    ("Wilkes", 7.1, 29),
    ("Wilkes", 7.1, 52),
    ("Recovery", 9.2, 48),
    ("Wilkes", 7.1, 12),
    ("Recovery", 9.2, 46),
    ("Aurora", 7.4, 10),
    ("Aurora", 7.4, 35),
    ("Pensacola-Pole", 15.0, 11),
    ("Gamburtsev-Transantarctic", 27.1, 28),
    ("Gamburtsev-Queen Maud", 28.1, 28),
]
x_text = sort(vcat(large_misi, small_misi, perimeter))
for i in eachindex(bifurcation_data)
    rectangle!(ax, (x_text[i] - 0.17, bifurcation_data[i][3] - 0.4), 0.4, bifurcation_data[i][2])
    text!(ax, x_text[i] + 0.2, bifurcation_data[i][3], text = bifurcation_data[i][1], rotation = π/2)
end

for k in 1:aqef.n_xps
    lines!(ax, aqef.f[k][1:s:end] ./ polar_amplification .+ f2020, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k], color = lcolor(cycling_colors[k]))
end

idx = sortperm(eql.f)
f = eql.f[idx]
V_sle = eql.V_sle[idx]
V_sle = filter_outliers(V_sle, mode = :high)
scatter!(ax, f ./ polar_amplification .+ f2020, V_sle;
    color = xpcolors["EQL"], label = "EQL", markersize = ms1)
axislegend(ax, position = :lb, nbanks = 1)

# Compute mean diff between eql and ref
aqef_itp = linear_interpolation(aqef.f[end] ./ polar_amplification .+ f2020, aqef.V_sle[end])
aqef_eql = aqef_itp.(f ./ polar_amplification .+ f2020)
mean_diff = mean(abs.(aqef_eql .- V_sle))
println("Mean difference REF vs EQL: $(round(mean_diff, digits=2)) m SLE")
lines(f, aqef_eql .- V_sle)

ax.xticks = 0:1:10
ax.xminorticks = 0:0.2:12
ax.yticks = 0:10:60
ax.yminorticks = IntervalsBetween(5)
ax.xminorgridvisible = true
ax.yminorgridvisible = true
ax.xlabel = L"GMT anomaly $f$ (K)"
ax.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(ax, 0, 60)
xlims!(ax, 1, 11)
fig1

axmap = Axis(fig1[1, 2], aspect = DataAspect())
cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)

ii = cropx+1:nx-cropx
jj = cropy+1:ny-cropy
XX = X[ii, jj]
YY = Y[ii, jj]

regions = [
    ("WAIS", (-1300, -300), 0, :black),
    ("EAIS", (900, -900), 0, :black),
]
subregions = [
    ("Antarctic \n Peninsula", (-2400, 1500), -π/5, :black),
    ("Amundsen", (-1600, -100), -π/3, :white),
    ("Siple Coast", (-700, -950), π/4, :white),
    ("Ross", (-200, -1100), 0, :white),
    ("Filchner- \n Ronne", (-1400, 400), 0, :white),
    ("Wilkes \n Basin", (700, -1950), 0, :black),
    ("Recovery \n Basin", (-450, 1000), 0, :black),
    ("Aurora \n Basin", (1700, -900), 0, :black),
    ("Vostok Lake", (900, -300), 0, :black),
    ("Gamburtsev", (500, 100), 0, :black),
    ("Amery \n Basin", (1000, 500), 0, :black),
    ("Transantarctic", (100, -300), -π/3, :black),
    ("Queen Maud Land", (-100, 1500), 0, :black),
    ("Kemp \n Land", (1500, 1000), 0, :black),
    ("Pensacola-Pole \n Basin", (-650, 400), -π/5, :black),
]
region_fontsize = 28
subregion_fontsize = 20


hidedecorations!(axmap)
heatmap!(axmap, xc, yc, z_bed; cmaps["z_bed2"]...)
heatmap!(axmap, xc, yc, z_srf .* f_ice, alpha = alpha; cmaps["z_srf"]...)
contour!(axmap, xc, yc, f_grnd .+ f_ice, levels = [1.9], color = :red, linewidth = 2)
for i in 1:3
    smooth_grline!(mask_ref)
end
contour!(axmap, xc, yc, mask_ref .== 2, levels = [0.5],
    color = :orange, linewidth = 2)
xlims!(axmap, extrema(XX))
ylims!(axmap, extrema(YY))
for i in eachindex(regions)
    text!(
        axmap,
        regions[i][2]...,
        text = regions[i][1],
        color = regions[i][4],
        fontsize = region_fontsize,
        font = :bold,
        rotation = regions[i][3],
    )
end
for i in eachindex(subregions)
    text!(
        axmap,
        subregions[i][2]...,
        text = subregions[i][1],
        color = subregions[i][4],
        fontsize = subregion_fontsize,
        font = :bold,
        rotation = subregions[i][3],
    )
end

relwidth = 0.35
Colorbar(fig1[2, 2], vertical = false, flipaxis = false, width = Relative(relwidth), halign = :left,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig1[2, 2], vertical = false, flipaxis = false, width = Relative(relwidth), halign = :right,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :orange, linewidth = 2)
Legend(fig1[2, 2], [elem_1, elem_2], ["Modelled", "Observed"], label = "Grounding zone", valign = 3)

text!(ax, 10.2, 56, color = :black, font = :bold, text="(a)", fontsize = 30)
text!(axmap, -2600, 2150, color = :black, font = :bold, text="(b)", fontsize = 30)
colsize!(fig1.layout, 1, 740)
colsize!(fig1.layout, 2, 800)
rowsize!(fig1.layout, 2, 50)
rowgap!(fig1.layout, 1, -70)
colgap!(fig1.layout, 1, -30)
fig1

save(plotsdir("v2/hysteresis/fig1-reduced.png"), fig1)
save(plotsdir("v2/hysteresis/fig1-reduced.pdf"), fig1)