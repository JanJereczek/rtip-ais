include("../../intro.jl")

T = Float32
visc_type = "lowvisc"    # "equil" or "aqef"
heatmap_frames = "aqef"    # "equil" or "aqef"

xps = [
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-normforcing-withrestarts"),
]
xp_labels = [
    "Quasi-equilibrium",
]
aqef = AQEFResults(T, xps)
eqldir = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
eql = EquilResults(T, eqldir)
lws = [5]
cycling_colors = [:steelblue1]

fn_ref = "/p/projects/megarun/ice_data/Antarctica/ANT-16KM/ANT-16KM_TOPO-RTOPO-2.0.1.nc"
mask_ref = ncread(fn_ref, "mask")
mask_ref[250:300, 150:200] .= 2

function smooth_grline!(mask)
    for I in CartesianIndices(mask)
        i, j = Tuple(I)
        if 1 < i < 380 && 1 < j < 380 &&
            mask[i+1, j] == 2 &&
            mask[i+1, j+1] == 2 &&
            mask[i+1, j-1] == 2 &&
            mask[i+1, j+1] == 2 &&
            mask[i+1, j-1] == 2 &&
            mask[i-1, j] == 2 &&
            mask[i-1, j+1] == 2 &&
            mask[i-1, j-1] == 2

            mask[i, j] = 2

        end
    end
end

for i in 1:3
    smooth_grline!(mask_ref)
end

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
xp_idx = aqef.n_xps
f2015 = 1.2

cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)
set_theme!(theme_latexfonts())
ms1, ms2 = 8, 18
nrows, ncols = 2, 1
fig = Figure(size=(750*0.9, 1050*0.9), fontsize = 24)
axs = [Axis(fig[i, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]
s = 400

shade = [(2.4, 2.5), (5.6, 5.95), (7.2, 7.4), (8.15, 8.3), (8.95, 9.1)]
# shade = [(1.2, 1.3), (4.6, 4.7), (6.1, 6.2), (7.0, 7.1), (7.8, 7.9)]
for i in eachindex(shade)
    if i == 1
        vlines!(axs[1, 1], (shade[i][1]:0.01:shade[i][2]), alpha = 0.2,
        color = :gray70, label = "Bifurcation", linewidth = 5)
    else
        vlines!(axs[1, 1], (shade[i][1]:0.01:shade[i][2]), alpha = 0.2, color = :gray70)
    end
end

text!(axs[1, 1], 11.2, 55, text = "a", color = :grey10, fontsize = 30, font = :bold)
axs[1, 1].xticks = 0:2:12
axs[1, 1].xminorticks = 0:0.2:12
axs[1, 1].yticks = 0:10:60
axs[1, 1].yminorticks = IntervalsBetween(10)
axs[1, 1].xminorgridvisible = true
axs[1, 1].yminorgridvisible = true
axs[1, 1].xaxisposition = :top
axs[1, 1].xlabel = L"GMT anomaly $f$ (K)"
axs[1, 1].ylabel = "AIS volume (m SLE)"
ylims!(axs[1, 1], 0, 60)
xlims!(axs[1, 1], 0, 12)
fig

file2D = joinpath(aqef.xps[xp_idx], "0", "yelmo2D.nc")
var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))
z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D, 1)

ii = cropx+1:nx-cropx
jj = cropy+1:ny-cropy
XX = X[ii, jj]
YY = Y[ii, jj]

region_labels = ["WAIS", "EAIS"]
region_positions = [(-1500, -500), (200, 0)]
subregion_labels = ["WSB", "RSB", "ASB"]
subregion_positions = [(600, -1750), (-400, 900), (1600, -800)]
region_color = :black
subregion_color = :honeydew
region_fontsize = 28
subregion_fontsize = 22

hidedecorations!(axs[2, 1])
heatmap!(axs[2, 1], xc, yc, z_bed; cmaps["z_bed2"]...)
heatmap!(axs[2, 1], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
contour!(axs[2, 1], xc, yc, f_grnd .+ f_ice, levels = [1.9],
    color = :red, linewidth = 2)
contour!(axs[2, 1], xc, yc, mask_ref .== 2, levels = [0.5],
    color = :orange, linewidth = 2)
xlims!(axs[2, 1], extrema(XX))
ylims!(axs[2, 1], extrema(YY))

for i in eachindex(region_labels)
    text!(axs[2, 1], region_positions[i]..., text = region_labels[i], color = region_color,
        fontsize = region_fontsize, font = :bold)
end
for i in eachindex(subregion_labels)
    text!(axs[2, 1], subregion_positions[i]..., text = subregion_labels[i], color = subregion_color,
        fontsize = subregion_fontsize, font = :bold)
end

for k in 1:aqef.n_xps
    lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification .+ f2015, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k], color = color = cycling_colors[k])
end

text!(axs[2, 1], 2300, 1800, text = "b", color = :white, fontsize = 30, font = :bold)

scatter!(axs[1, 1], eql.f ./ polar_amplification .+ f2015, eql.V_sle;
    color = :black, label = "Equilibrium", markersize = ms1)

axislegend(axs[1, 1], position = :lb, nbanks = 1)
relheight = 0.8
Colorbar(fig[2, 0], vertical = true, height = Relative(relheight), flipaxis = false,
    label = "Bed elevation (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[2, 2], vertical = true, height = Relative(relheight),
    label = "Ice surface elevation (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
l1 = LineElement(color = :orange, linewidth = 2)
l2 = LineElement(color = :red, linewidth = 2)
Legend(fig[3, :], [l1, l2], ["Observed", "Modelled grounding line"],
    nbanks = 2, valign = :bottom)
    
rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
colgap!(fig.layout, 1, -60)
rowsize!(fig.layout, 3, 40)

text!(axs[1, 1], 2.4, 20, text = "WAIS", rotation = π/2)
text!(axs[1, 1], 5.6, 15, text = "low lat. WSB", rotation = π/2)
text!(axs[1, 1], 7.2, 20, text = "RSB", rotation = π/2)
text!(axs[1, 1], 8.15, 35, text = "high lat. WSB", rotation = π/2)
text!(axs[1, 1], 8.95, 40, text = "ASB", rotation = π/2)
save(plotsdir("16km/hysteresis/retreat-minimal-$visc_type.png"), fig)
save(plotsdir("16km/hysteresis/retreat-minimal-$visc_type.pdf"), fig)