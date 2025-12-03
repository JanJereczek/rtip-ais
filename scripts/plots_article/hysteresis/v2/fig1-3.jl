include("../../../intro.jl")

T = Float32
heatmap_frames = "aqef"    # "equil" or "aqef"
xps = [
    datadir("output/ais/v2/hyster/retreat/aqef/minvisc/ocn"),
    datadir("output/ais/v2/hyster/retreat/aqef/minvisc/atm"),
    datadir("output/ais/v2/hyster/retreat/aqef/minvisc/dpr"),
    datadir("output/ais/v2/hyster/retreat/aqef/minvisc/refnomslow"),
]
xp_labels = [
    "OCN",
    "ATM",
    "DPR",
    "REF",
]
aqef = AQEFResults(T, xps)
eqldir = datadir("output/ais/v2/hyster/retreat/equil/refnomslow")
eql = EquilResults(T, eqldir)
lws = [5, 5, 5, 10]
cycling_colors = [xpcolors[l] for l in xp_labels]
cycling_colors[end] = :gray10

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
xp_idx = aqef.n_xps
f2015 = 1.2

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
perimeter = [9.7, 10.1]
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
    ("Gamburtsev-Dronning Maud", 28.1, 28),
]
x_text = sort(vcat(large_misi, small_misi, perimeter))
for i in eachindex(bifurcation_data)
    rectangle!(ax, (x_text[i] - 0.17, bifurcation_data[i][3] - 0.4), 0.4, bifurcation_data[i][2])
    text!(ax, x_text[i] + 0.2, bifurcation_data[i][3], text = bifurcation_data[i][1], rotation = π/2)
end

for k in 1:aqef.n_xps
    lines!(ax, aqef.f[k][1:s:end] ./ polar_amplification .+ f2015, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k], color = lcolor(cycling_colors[k]))
end

scatter!(ax, eql.f ./ polar_amplification .+ f2015, eql.V_sle;
    color = xpcolors["REF"], label = "EQL", markersize = ms1)

axislegend(ax, position = :lb, nbanks = 1)

ax.xticks = 0:2:12
ax.xminorticks = 0:0.2:12
ax.yticks = 0:10:60
ax.yminorticks = IntervalsBetween(10)
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
    ("EAIS", (900, -500), 0, :black),
    ("APIS", (-2200, 1300), -π/5, :black),
]
subregions = [
    ("Amundsen", (-1600, -100), -π/3, :white),
    ("Ross", (-200, -1100), 0, :white),
    ("Filchner- \n Ronne", (-1400, 400), 0, :white),
    ("Wilkes", (700, -1750), 0, :black),
    ("Recovery", (-450, 1100), 0, :black),
    ("Aurora", (1600, -1000), 0, :black),
    ("Gamburtsev", (500, 300), 0, :black),
    ("Amery", (1400, 600), 0, :white),
    ("Transantarctic", (100, -300), -π/3, :black),
    ("Dronning Maud", (0, 1500), 0, :black),
    ("Pensacola-Pole", (-600, 400), -π/5, :black),
]
region_fontsize = 30
subregion_fontsize = 22


hidedecorations!(axmap)
heatmap!(axmap, xc, yc, z_bed; cmaps["z_bed2"]...)
heatmap!(axmap, xc, yc, z_srf .* f_ice, alpha = alpha; cmaps["z_srf"]...)
contour!(axmap, xc, yc, f_grnd .+ f_ice, levels = [1.9], color = :red, linewidth = 2)
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
fig1
colsize!(fig1.layout, 1, 740)
colsize!(fig1.layout, 2, 800)
rowsize!(fig1.layout, 2, 50)
rowgap!(fig1.layout, 1, -70)
colgap!(fig1.layout, 1, -30)
save(plotsdir("v2/hysteresis/fig1.png"), fig1)
save(plotsdir("v2/hysteresis/fig1.pdf"), fig1)

################################
# Fig 2
################################
nrows, ncols = 2, 4
fframes = sort(vcat(large_misi, large_misi .+ 0.1))
forcing_frames = transpose(reshape(fframes, ncols, nrows))
state_labels = ["a" "b" "c" "d"; "e" "f" "g" "h"]
fig2 = Figure(size=(1600, 850), fontsize = 24)
axs = [Axis(fig2[i+1, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]

xl = [
    (-1900, -900),
    (600, 1500),
    (-600, 400),
    (1400, 2600)
]
yl = [
    (-800, 200),
    (-2300, -1400),
    (800, 1800),
    (-1300, -100),
]

xlims_frames = permutedims(reshape(repeat(xl, inner=2), ncols, nrows))
ylims_frames = permutedims(reshape(repeat(yl, inner=2), ncols, nrows))

for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    forcing = forcing_frames[i, j]

    if heatmap_frames == "aqef"
        i3 = findfirst(aqef.f[xp_idx] ./ polar_amplification .+ f2015 .>= forcing)
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
    @show i3, f_eq, V_eq, frame_index

    hidedecorations!(axs[i, j])
    heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
    heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
    contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
        color = :red, linewidth = 2)
    contour!(axs[i, j], xc, yc, (xlims_frames[i, j][1] .< X .< xlims_frames[i, j][2]) .&
        (ylims_frames[i, j][1] .< Y .< ylims_frames[i, j][2]), levels = [0.5],
        color = :darkred, linewidth = 3)
    statlab = state_labels[i, j]
    text!(axs[i, j], -2500, -2450, color = :black, font = :bold,
        text="("*statlab*") f = $(round(forcing, digits=2)) K", fontsize = 24)
    xlims!(axs[i, j], extrema(XX))
    ylims!(axs[i, j], extrema(YY))
end

relwidth = 0.4
Colorbar(fig2[1, 1:2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig2[1, 3:4], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :darkred, linewidth = 2)
Legend(fig2[1, 2:3], [elem_1, elem_2], ["Grounding line", "Highlighted region"], valign = :bottom)

rowsize_base = 350
rowgap!(fig2.layout, 5)
colgap!(fig2.layout, 5)
rowgap!(fig2.layout, 1, 30)
rowsize!(fig2.layout, 1, 1)
rowsize!(fig2.layout, 2, rowsize_base)
rowsize!(fig2.layout, 3, rowsize_base)
colsize!(fig2.layout, 1, rowsize_base*aratio)
colsize!(fig2.layout, 2, rowsize_base*aratio)
colsize!(fig2.layout, 3, rowsize_base*aratio)
colsize!(fig2.layout, 4, rowsize_base*aratio)
fig2
save(plotsdir("v2/hysteresis/fig2.png"), fig2)
save(plotsdir("v2/hysteresis/fig2.pdf"), fig2)


################################
# Fig 3
################################

xp = joinpath(aqef.xps[end], "0")
file1D = joinpath(xp, "yelmo1D.nc")
file2Dsm = joinpath(xp, "yelmo2Dsm.nc")
file2D = joinpath(xp, "yelmo2D.nc")

t1D = ncread(file1D, "time")
t2Dsm = ncread(file2Dsm, "time")
t2D = ncread(file2D, "time")
dt_2D = mean(diff(t2D))
f = aqef.f[end] ./ polar_amplification .+ f2015

n_grz = 7       # number of plotted grounding zones
di_grz = [2, 4, 6, 4]
dt_grz = di_grz .* dt_2D
n_tra = length(di_grz)
f_bif = large_misi
i_bif = [findlast(f .<= f_bif[i]) for i in 1:n_tra] .+ 2
t_bif = [t1D[i_bif[i]] for i in 1:n_tra]
t_end = t_bif .+ dt_grz .* n_grz

mask(X, Y, xl, yl) = (xl[1] .<= X .<= xl[2]) .& (yl[1] .<= Y .<= yl[2])
masks = [ mask(X, Y, xl[i], yl[i]) for i in 1:n_tra ]
mbs = [ masked_massbalance(file2Dsm, t_bif[i]- dt_2D, t_end[i], masks[i]) for i in 1:n_tra ]

img = FileIO.load(datadir("processed/transects-v2.png"))
rgba2gray_smooth(x) = x.r + x.g + x.b
img = [rgba2gray_smooth(img[i, j]) for i in 1:size(img, 1), j in 1:size(img, 2)]
mask_vals = sort(unique(img))
mask_order = [2, 3, 4, 5]
mask_transects = rotr90(img)
transect_begin = [:left, :bottom, :left, :right]
region_labels = ["Amundsen", "Wikes", "Recovery", "Aurora"]

fs = 20
lw1 = 3
lw2 = 5
fig3 = Figure(size = (1050, 1200), fontsize = fs)
axs_ts = [Axis(fig3[i, 1], aspect = AxisAspect(0.85)) for i in 2:n_tra+1]
axs_hm = [Axis(fig3[i, 2], aspect = AxisAspect(1)) for i in 2:n_tra+1]
axs_tr = [Axis(fig3[i, 3:4], aspect = AxisAspect(2)) for i in 2:n_tra+1]
catjet = cgrad(:jet, range(0, stop = 1, length = n_grz+1), categorical = true)
tsjet = cgrad(:jet)
extract_idx(idx) = Tuple(idx.I)
tr_color = :black   # :magenta

for i in eachindex(i_bif)

    axs_ts[i].ylabel = "Mean MB (m/yr)"

    hideydecorations!(axs_hm[i])
    axs_hm[i].xticksvisible = false
    axs_hm[i].yticksvisible = false
    axs_hm[i].xticklabelsvisible = false
    axs_hm[i].yticklabelsvisible = false
    axs_hm[i].xlabel = region_labels[i]

    axs_tr[i].xgridvisible = false
    axs_tr[i].ygridvisible = false
    axs_tr[i].yaxisposition = :right
    axs_tr[i].ylabel = "Elevation (km)"
    axs_tr[i].yticks = latexifyticks(-2:1:3, 1e3)
    axs_tr[i].bottomspinecolor = tr_color
    axs_tr[i].topspinecolor = tr_color
    axs_tr[i].spinewidth = 3
    
    t = t1D[i_bif[i]]
    k = argmin(abs.(t .- t2D))
    heatmap!(axs_hm[i], xc, yc, ncslice(file2D, "z_bed", k); cmaps["z_bed6"]...)
        grz_steps = range(k, step = di_grz[i], length = n_grz)
    vlines!(axs_ts[i], t2D[grz_steps] ./ 1f3, color = [(c, 0.3) for c in catjet],
        linewidth = 6)
    hlines!(axs_ts[i], 0, color = :black, linewidth = lw1, linestyle = :dash)
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].bmb, linewidth = lw2,
        label = "basal")
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].smb, linewidth = lw2,
        label = "surface")
    # lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].cmb, linewidth = lw2,
    #     label = "calving", color = :gray70)
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].mb_net,
        color = :gray30,
        label = "net",
        # color = mbs[i].time, colormap = tsjet,
        linewidth = lw2)

    tr_mask = mask_transects .== mask_vals[mask_order[i]]
    ordered_idx, rt = indices_from_mask(tr_mask, X, Y, transect_begin[i])
    ii = [extract_idx(ordered_idx[i])[1] for i in 1:length(ordered_idx)]
    jj = [extract_idx(ordered_idx[i])[2] for i in 1:length(ordered_idx)]
    plot_simple_transect(axs_tr[i],
        ncslice(file2D, "z_bed", k),
        ncslice(file2D, "f_grnd", k),
        ncslice(file2D, "H_ice", k),
        ncslice(file2D, "z_sl", k),
        ii[1:end-1],
        jj[1:end-1],
        rt[1:end-1],
    )

    for l in eachindex(grz_steps)
        contour!(axs_hm[i], xc, yc,
            ncslice(file2D, "f_grnd", grz_steps[l]) .*
            (ncslice(file2D, "H_ice", grz_steps[l]) .> 100),
            levels = [0.5], color = catjet[l], linewidth = lw1)
        z_bed = ncslice(file2D, "z_bed", grz_steps[l])
        H_ice = ncslice(file2D, "H_ice", grz_steps[l])
        f_grnd = ncslice(file2D, "f_grnd", grz_steps[l])
        z_srf = z_bed .+ H_ice .* f_grnd
        no_ice = (f_grnd .< 0.01) .|| (H_ice .< 1)
        z_srf[no_ice] .= NaN
        z_srf_cross = crossection(z_srf, ii[1:end-1], jj[1:end-1])
        z_bed_cross = crossection(z_bed, ii[1:end-1], jj[1:end-1])
        grindex = findlast(isnan.(z_srf_cross))
        z_srf_cross[grindex] = z_bed_cross[grindex]

        # if (length(z_srf_cross[.!isnan.(z_srf_cross)])) > 1 &&
        #     (z_srf_cross[.!isnan.(z_srf_cross)][1] > z_srf_cross[.!isnan.(z_srf_cross)][2])
            
        #     z_srf_cross[.!isnan.(z_srf_cross)][1] = NaN
        # end

        if i == 1 && l == 4
            z_srf_cross[grindex-1] = NaN
        end
        
        lines!(axs_tr[i], rt[1:end-1], z_srf_cross, color = catjet[l], linewidth = lw1)
    end

    heatmap!(axs_hm[i], xc, yc, tr_mask, colormap = cgrad([tr_color, tr_color]),
        colorrange = (1, 1.1), lowclip = :transparent, highclip = :transparent)
    xlims!(axs_hm[i], xl[i])
    ylims!(axs_hm[i], yl[i])
end

axs_ts[end].xlabel = "Time (kyr)"
axs_tr[end].xlabel = "Distance along transect (km)"
Legend(fig3[1, 1], axs_ts[end], nbanks = 2, fontsize = fs-4, valign = :bottom)
Colorbar(fig3[1, 2], label = "Bed elevation (km)", vertical = false,
    width = Relative(0.9), ticks = latexifyticks(-1:1, 1e3), flipaxis = true;
    cmaps["z_bed6"]...)
Legend(fig3[1, 3:4], axs_tr[end], nbanks = 2, fontsize = fs-4, valign = :bottom)
# le = LineElement(
#     points = Point2f.(range(0, 1, 10), 0.5 .+ 0.5 * sin.(range(0, 2pi, 10))),
#     color = 1:10,
#     colormap = tsjet,
# )
# Legend(fig3[1, 1], [le], ["net"], framevisible = false, halign = :right, valign = :bottom, patchlabelgap = 10)
rowsize!(fig3.layout, 1, 10)
colgap!(fig3.layout, 5)
rowgap!(fig3.layout, -15)
rowgap!(fig3.layout, 1, -5)
fig3

save(plotsdir("v2/hysteresis/fig3.png"), fig3)
save(plotsdir("v2/hysteresis/fig3.pdf"), fig3)