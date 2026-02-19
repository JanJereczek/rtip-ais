include("../../../intro.jl")

T = Float32
dir = datadir("output/ais/v2/hyster/retreat/aqef/minvisc")
xps = [
    "$dir/refm2slow",
    "$dir/refm1slow",
    "$dir/refp1slow",
    "$dir/refp2slow",
    "$dir/refp2slow-restarted",
    "$dir/refnomslow",
]

xp_labels = [
    L"$-2 \, \sigma$",
    L"$-1 \, \sigma$",
    L"$+1 \, \sigma$",
    L"$+2 \, \sigma$",
    nothing,
    L"$0 \, \sigma$",
]
aqef = AQEFResults(T, xps)

# eqldir1 = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
# eql1 = EquilResults(T, eqldir1)
# eqldir2 = datadir("output/ais/hyster/16km/regrowth/equil")
# eql2 = EquilResults(T, eqldir2)

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

lw1, lw2, lw3 = 2, 3, 3
ms1, ms2 = 15, 20
lws = vcat(fill(lw2, aqef.n_xps - 1), [lw3])
cycling_colors = [
    viscmap[1],
    viscmap[2],
    viscmap[4],
    viscmap[5],
    viscmap[5],
    viscmap[3],
]

s = 50
uniqueidx(v) = unique(i -> v[i], eachindex(v))
file2D = joinpath(aqef.xps[4], "0", "yelmo2D.nc")
var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))
z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D, 1)

##################################
# Fig
##################################

set_theme!(theme_latexfonts())
cropx, cropy = 20, 38
aratio = (381 - 2*cropx) / (381 - 2*cropy)
fig = Figure(size=(700, 1050), fontsize = 24)
ax1 = Axis(fig[1, 1], aspect = AxisAspect(aratio))
# ax2 = Axis(fig[3, 1])
ax3 = Axis(fig[2, 1], aspect = DataAspect()) #AxisAspect(aratio)

####################################
# Ax1
#####################################
f_gmt_bif = [1.98, 5.87, 7.1, 8.0, 8.65]
shade = [(1.9, 2.0), (5.85, 5.95), (6.8, 7.1), (7.85, 7.95), (8.6, 8.7)]
for i in eachindex(shade)
    vlines!(ax1, (shade[i][1]:0.01:shade[i][2]), alpha = 0.2, color = :gray70)
end
vlines!(ax1, 1f6, color = :gray70, label = "Bifurcation", linewidth = 5) 
text!(ax1, 1.9, 40, text = "WAIS", rotation = π/2)
text!(ax1, 5.8, 25, text = "WSB", rotation = π/2)
text!(ax1, 6.8, 20, text = "RSB", rotation = π/2)
text!(ax1, 7.8, 6, text = "ASB ext.", rotation = π/2)
text!(ax1, 9.3, 20, text = "ASB int.", rotation = π/2)

V_bif = [54, 42, 31, 20, 13]
for k in 1:aqef.n_xps
    idx = uniqueidx(aqef.t_1D[k])
    lines!(
        ax1,
        aqef.f[k][idx][1:s:end] ./ polar_amplification .+ f2020,
        aqef.V_sle[k][idx][1:s:end],
        linewidth = lws[k],
        label = xp_labels[k],
        color = cycling_colors[k],
    )
    i_bif = [findlast(aqef.V_sle[k][idx] .> V_bif[i]) for i in eachindex(V_bif)]
    i_bif[i_bif .== nothing] .= length(aqef.V_sle[k][idx])
    f_rtip = [aqef.f[k][idx][i_bif[i]] for i in eachindex(i_bif)]
    @show f_rtip
end

text!(ax1, 11.2, 55, text = "a", color = :grey10, fontsize = 30, font = :bold)
ax1.xticks = 0:2:12
ax1.xminorticks = 0:0.2:12
ax1.yticks = 0:10:60
ax1.yminorticks = 0:2:60
ax1.xminorgridvisible = true
ax1.yminorgridvisible = true
ax1.xaxisposition = :top
ax1.xlabel = L"GMT anomaly $f$ (K)"
ax1.ylabel = "AIS volume (m SLE)"
ylims!(ax1, 0, 60)
xlims!(ax1, 0, 12)
axislegend(ax1, position = :lb)

################################
# Ax3
#################################
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

hidedecorations!(ax3)
heatmap!(ax3, xc, yc, z_bed; cmaps["z_bed2"]...)
heatmap!(ax3, xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
contour!(ax3, xc, yc, f_grnd .+ f_ice, levels = [1.9],
    color = :red, linewidth = lw1)
contour!(ax3, xc, yc, mask_ref .== 2, levels = [0.5],
    color = :orange, linewidth = lw1)
xlims!(ax3, extrema(XX))
ylims!(ax3, extrema(YY))

for i in eachindex(region_labels)
    text!(ax3, region_positions[i]..., text = region_labels[i], color = region_color,
        fontsize = region_fontsize, font = :bold)
end
for i in eachindex(subregion_labels)
    text!(ax3, subregion_positions[i]..., text = subregion_labels[i], color = subregion_color,
        fontsize = subregion_fontsize, font = :bold)
end

relheight = 0.8
Colorbar(fig[2, 0], vertical = true, height = Relative(relheight), flipaxis = false,
    label = "Bed elevation (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[2, 2], vertical = true, height = Relative(relheight),
    label = "Ice surface elevation (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
l1 = LineElement(color = :orange, linewidth = lw1)
l2 = LineElement(color = :red, linewidth = lw1)
Legend(fig[3, :], [l1, l2], ["Observed", "Modelled grounding line"],
    nbanks = 2, valign = :bottom)

text!(ax1, 11.2, 55, text = "a", color = :grey10, fontsize = 30, font = :bold)
text!(ax3, 2300, 1800, text = "b", color = :white, fontsize = 30, font = :bold)

rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
colgap!(fig.layout, 1, -60)
rowsize!(fig.layout, 3, 40)
fig

save(plotsdir("v2/rtip/btip.png"), fig)