include("../intro.jl")

T = Float32
polar_amplification = 1.8
f_to = 0.25

xp = datadir("output/ais/hyster/16km/retreat/aqef/pmpt-lowvisc-normforcing-withrestarts/0")
file1D = joinpath(xp, "yelmo1D.nc")
file2Dsm = joinpath(xp, "yelmo2Dsm.nc")
file2D = joinpath(xp, "yelmo2D.nc")

t1D = ncread(file1D, "time")
t2Dsm = ncread(file2Dsm, "time")
t2D = ncread(file2D, "time")
f = ncread(file1D, "hyst_f_now") ./ polar_amplification
x = ncread(file2Dsm, "xc")
y = ncread(file2Dsm, "yc")
X, Y = ndgrid(x, y)

n_grz = 7       # number of plotted grounding zones
n_tra = 5       # number of plotted transects
di_grz = [1, 3, 5, 7, 3]
dt_2D = mean(diff(t2D))
dt_grz = di_grz .* dt_2D
f_bif = [1.2, 4.5, 6.0, 6.9, 7.8]
i_bif = [findlast(f .<= f_bif[i]) for i in 1:n_tra] .+ 2
t_bif = [t1D[i_bif[i]] for i in 1:n_tra]
t_end = t_bif .+ dt_grz .* n_grz

# xl = [
#     (-2000, -1000),
#     (700, 1400),
#     (-500, 500),
#     (400, 1400),
#     (1300, 2500)
# ]
# yl = [
#     (-600, 400),
#     (-2300, -1600),
#     (900, 1900),
#     (-1700, -700),
#     (-1100, 100),
# ]

xl = [
    (-1900, -900),
    (700, 1400),
    (-600, 400),
    (400, 1400),
    (1300, 2500)
]
yl = [
    (-800, 200),
    (-2300, -1600),
    (800, 1800),
    (-1600, -600),
    (-1000, 0),
]
mask(X, Y, xl, yl) = (xl[1] .<= X .<= xl[2]) .& (yl[1] .<= Y .<= yl[2])
masks = [ mask(X, Y, xl[i], yl[i]) for i in 1:n_tra ]
mbs = [ masked_massbalance(file2Dsm, t_bif[i]- dt_2D, t_end[i], masks[i]) for i in 1:n_tra ]

img = FileIO.load(datadir("processed/transects.png"))
rgba2gray_smooth(x) = x.r + x.g + x.b
img = [rgba2gray_smooth(img[i, j]) for i in 1:size(img, 1), j in 1:size(img, 2)]
mask_vals = sort(unique(img))
mask_order = [2, 4, 5, 3, 6]
mask_transects = rotr90(img)
transect_begin = [:left, :bottom, :left, :bottom, :right]
region_labels = ["ASE", "WSB 1", "RSB", "WSB 2", "ASB"]

set_theme!(theme_latexfonts())
fs = 20
lw1 = 3
lw2 = 5
fig = Figure(size = (2100 / 2, 2900 / 2), fontsize = fs)
axs_ts = [Axis(fig[i, 1], aspect = AxisAspect(0.85)) for i in 2:n_tra+1]
axs_hm = [Axis(fig[i, 2], aspect = AxisAspect(1)) for i in 2:n_tra+1]
axs_tr = [Axis(fig[i, 3:4], aspect = AxisAspect(2)) for i in 2:n_tra+1]
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
    heatmap!(axs_hm[i], x, y, ncslice(file2D, "z_bed", k); cmaps["z_bed6"]...)
        grz_steps = range(k, step = di_grz[i], length = n_grz)
    vlines!(axs_ts[i], t2D[grz_steps] ./ 1f3, color = [(c, 0.8) for c in catjet],
        linewidth = 1)
    hlines!(axs_ts[i], 0, color = :black, linewidth = lw1, linestyle = :dash)
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].bmb, linewidth = lw1,
        label = "basal")
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].smb, linewidth = lw1,
        label = "surface")
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].cmb, linewidth = lw1,
        label = "calving", color = :gray70)
    lines!(axs_ts[i], mbs[i].time ./ 1f3, mbs[i].mb_net,
        color = :gray30,
        label = "net",
        # color = mbs[i].time, colormap = tsjet,
        linewidth = lw1)

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
        contour!(axs_hm[i], x, y,
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

    heatmap!(axs_hm[i], x, y, tr_mask, colormap = cgrad([tr_color, tr_color]),
        colorrange = (1, 1.1), lowclip = :transparent, highclip = :transparent)
    xlims!(axs_hm[i], xl[i])
    ylims!(axs_hm[i], yl[i])
end

axs_ts[end].xlabel = "Time (kyr)"
axs_tr[end].xlabel = "Distance along transect (km)"
Legend(fig[1, 1], axs_ts[end], nbanks = 2, fontsize = fs-4, valign = :bottom)
Colorbar(fig[1, 2], label = "Bed elevation (km)", vertical = false,
    width = Relative(0.9), ticks = latexifyticks(-1:1, 1e3), flipaxis = true;
    cmaps["z_bed6"]...)
Legend(fig[1, 3:4], axs_tr[end], nbanks = 2, fontsize = fs-4, valign = :bottom)
# le = LineElement(
#     points = Point2f.(range(0, 1, 10), 0.5 .+ 0.5 * sin.(range(0, 2pi, 10))),
#     color = 1:10,
#     colormap = tsjet,
# )
# Legend(fig[1, 1], [le], ["net"], framevisible = false, halign = :right, valign = :bottom, patchlabelgap = 10)
rowsize!(fig.layout, 1, 10)
colgap!(fig.layout, 5)
rowgap!(fig.layout, -15)
rowgap!(fig.layout, 1, -5)
save(plotsdir("16km/hysteresis/transects.png"), fig)