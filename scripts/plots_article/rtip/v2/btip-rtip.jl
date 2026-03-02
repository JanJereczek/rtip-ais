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

eqldir = datadir("output/ais/v2/hyster/retreat/equil/refnomslow")
eql = EquilResults(T, eqldir)

sr = StepRtip{Float32}(visc_cases)
prefix = "output/ais/v2/ramps/"
dirs = [
    datadir("$prefix/steps-sigmarange"),
]
for dir in dirs
    aggregate_xp!(sr, dir)
end
sort!(sr)

fn_ref = "/p/projects/megarun/ice_data/Antarctica/ANT-16KM/ANT-16KM_TOPO-RTOPO-2.0.1.nc"
mask_ref = ncread(fn_ref, "mask")
mask_ref[250:300, 150:200] .= 2
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
file2D = joinpath(aqef.xps[4], "0", "yelmo2D.nc")
var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = length(xc), length(yc)
z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D, 1)

##################################
# Fig
##################################

set_theme!(theme_latexfonts())
cropx, cropy = 20, 38
aratio = (nx - 2*cropx) / (ny - 2*cropy)
fig = Figure(size = (1600, 600), fontsize = 24)
ax_ssp = Axis(fig[1, 1], aspect = AxisAspect(aratio))
ax1 = Axis(fig[1, 1], aspect = AxisAspect(aratio))
ax2 = Axis(fig[1, 2], aspect = DataAspect())
ax3 = Axis(fig[1, 3], aspect = AxisAspect(aratio))
ax4 = Axis(fig[1, 3], aspect = AxisAspect(aratio))
rw = 0.4

hidedecorations!(ax3)
hidedecorations!(ax_ssp)
ylims!(ax1, 0, 60)
xlims!(ax1, 1, 11)
xlims!(ax_ssp, 1, 11)
vlines!(ax1, 8, color = :black, linestyle = :dash, label = "RCP8.5, 2300", linewidth = 3)

####################################
# Ax1
#####################################

V_bif = [54, 42, 31, 20, 13]
f_bif_atm = []
for k in 1:aqef.n_xps
    idx = uniqueidx(aqef.t_1D[k])
    i_bif = [findlast(aqef.V_sle[k][idx] .> V_bif[i]) for i in eachindex(V_bif)]
    f_bif_elem = [missingselect(aqef.f[k][idx], i_bif[i]) for i in eachindex(i_bif)]
    push!(f_bif_atm, f_bif_elem)
end
f_bif_atm = [f_bif_atm[1], f_bif_atm[2], f_bif_atm[6], f_bif_atm[3], skipmissingmin.(f_bif_atm[4], f_bif_atm[5])]
f_bif = [round.(f_bif_atm[i] ./ polar_amplification .+ f2020, digits = 2) for i in eachindex(f_bif_atm)]
@show f_bif

f_bif_wais = [f_bif[i][1] for i in eachindex(f_bif)]
f_bif_wsb = [f_bif[i][2] for i in eachindex(f_bif)]
f_bif_rsb = [f_bif[i][3] for i in eachindex(f_bif)]
f_bif_asb_ext = [f_bif[i][4] for i in eachindex(f_bif)]
f_bif_asb_int = [f_bif[i][5] for i in eachindex(f_bif)]

shade = [extrema(skipmissing(f)) for f in [f_bif_wais, f_bif_wsb, f_bif_rsb, f_bif_asb_ext, f_bif_asb_int]]
for i in eachindex(shade)
    vlines!(ax1, (shade[i][1]:0.1:shade[i][2]), alpha = 0.5, color = :mediumpurple, linewidth = 5)
end
vlines!(ax1, 1f6, color = :mediumpurple, label = "Bifurcation", linewidth = 5)

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
end

idx = sortperm(eql.f)
f = eql.f[idx]
V_sle = eql.V_sle[idx]
V_sle = filter_outliers(V_sle, mode = :high)
scatter!(
    ax1,
    f ./ polar_amplification .+ f2020,
    V_sle,
    color = :black,
    markersize = 6,
    label = "Equilibrium",
)

ssp2100 = [ssp126_2100, ssp245_2100, ssp370_2100, ssp585_2100]
ssp2100_hi = [ssp126_2100_hi, ssp245_2100_hi,
    ssp370_2100_hi, ssp585_2100_hi]
ssp_labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
ssp_short_labels = [
    L"\bar{f}_\text{1-2.6}",
    L"\bar{f}_\text{2-4.5}",
    L"\bar{f}_\text{3-7.0}",
    L"\bar{f}_\text{5-8.5}",
]
ssphi_short_labels = [
    L"\hat{f}_\text{1-2.6}",
    L"\hat{f}_\text{2-4.5}",
    L"\hat{f}_\text{3-7.0}",
    L"\hat{f}_\text{5-8.5}",
]
ssp_line_opts = (linewidth = lw1, alpha = 0.7)
ssphi_line_opts = (linestyle = :dash, linewidth = 3, alpha = 0.7)
# for (i, ssp) in enumerate(ssp2100)
#     vlines!(ax_ssp, ssp, label = ssp_short_labels[i], color = xpcolors[ssp_labels[i]];
#         ssp_line_opts...)
# end
# for (i, ssphi) in enumerate(ssp2100_hi)
#     vlines!(ax_ssp, ssphi, label = ssphi_short_labels[i], color = xpcolors[ssp_labels[i]];
#         ssphi_line_opts...)
# end

bifurcation_data = [
    ("WAIS", 11, 40),
    ("Wilkes", 12, 20),
    ("Recovery", 15, 40),
    ("Aurora ext.", 18.5, 30),
    ("Aurora int.", 18.5, 20),
]
x_text = f_bif[end]
for i in eachindex(bifurcation_data)
    rectangle!(ax1, (x_text[i] - 0.3, bifurcation_data[i][3] - 0.6), 0.65, bifurcation_data[i][2])
    text!(ax1, x_text[i] + 0.3, bifurcation_data[i][3], text = bifurcation_data[i][1], rotation = π/2)
end

axislegend(ax1, position = :lb)
# Legend(fig[0, 1], ax_ssp, nbanks = 4, valign = :bottom)
text!(ax1, 10.2, 53, text = "a", color = :grey10, fontsize = 30, font = :bold)
ax1.xticks = 0:1:10
ax1.xminorticks = 0:0.2:12
ax1.yticks = 0:10:60
ax1.yminorticks = 0:2:60
ax1.xaxisposition = :bottom
ax1.xlabel = "GMT anomaly (K)"
ax1.ylabel = L"AIS volume ($\mathrm{m \, SLE}$)"

################################
# Ax2
#################################
ii = cropx+1:nx-cropx
jj = cropy+1:ny-cropy
x = view(xc, ii)
y = view(yc, jj)
region_labels = ["WAIS", "EAIS"]
region_positions = [(-1500, -500), (200, 0)]
subregion_labels = ["Wilkes", "Recovery", "Aurora"]
subregion_positions = [(600, -1750), (-400, 900), (1400, -800)]
region_color = :black
subregion_color = :black
region_fontsize = 28
subregion_fontsize = 22

hidedecorations!(ax2)
heatmap!(ax2, xc, yc, z_bed; cmaps["z_bed2"]...)
heatmap!(ax2, xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
contour!(ax2, xc, yc, f_grnd .+ f_ice, levels = [1.9],
    color = :red, linewidth = lw1)
contour!(ax2, xc, yc, mask_ref .== 2, levels = [0.5],
    color = :orange, linewidth = lw1)
xlims!(ax2, extrema(x))
ylims!(ax2, extrema(y))

for i in eachindex(region_labels)
    text!(ax2, region_positions[i]..., text = region_labels[i], color = region_color,
        fontsize = region_fontsize, font = :bold)
end
for i in eachindex(subregion_labels)
    text!(ax2, subregion_positions[i]..., text = subregion_labels[i], color = subregion_color,
        fontsize = subregion_fontsize, font = :bold)
end

Colorbar(
    fig[0, 2],
    vertical = false,
    width = Relative(rw),
    flipaxis = true,
    label = "Bed elevation (km)",
    halign = :left,
    ticks = latexifyticks(-6:2:2, 1f3);
    cmaps["z_bed2"]...,
)
Colorbar(
    fig[0, 2],
    vertical = false,
    flipaxis = true,
    width = Relative(rw),
    label = "Ice surface elevation (km)",
    halign = :right,
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4));
    cmaps["z_srf"]...,
)
l1 = LineElement(color = :orange, linewidth = lw1)
l2 = LineElement(color = :red, linewidth = lw1)
Legend(fig[2, 2], [l1, l2], ["Observed", "Modelled grounding line"],
    nbanks = 2, valign = :bottom)

text!(ax2, 2300, 1800, text = "b", color = :white, fontsize = 30, font = :bold)

####################################
# Ax4
####################################

basins = ["WAIS", "Wilkes", "Recovery", "Aurora\n ext.", "Aurora\n int."]
i_rtip = [[findlast(V .> V_bif[i]) + 1 for i in eachindex(V_bif)] for V in sr.V]
f_rtip = [sr.f[i][i_rtip[i]] for i in eachindex(i_rtip)]

# Based on the full WSB R-tip diagram, the isolated dots of tipping can be attributed to chaos + noise and should be ignored.
for i in eachindex(f_rtip)
    f_rtip[i][2] = max.(f_rtip[i][2], 5.755)
end

f_bif_cat = vcat(f_bif...)
f_bif_cat[ismissing.(f_bif_cat)] .= NaN
r_gap = vcat([f_rtip[i] .- f_bif[i] for i in eachindex(f_rtip)]...)
r_gap = round.(r_gap, digits = 2)
b_gap = vcat([f_rtip[i] .- f_bif[1] for i in eachindex(f_rtip)]...)
b_gap = round.(b_gap, digits = 2)

tbl_gap = (
    cat = repeat(1:5, outer = 5),
    height = r_gap,
    grp = repeat(1:5, inner = 5),
)
tbl_bgap = (
    cat = repeat(1:5, outer = 5),
    height = b_gap,
    grp = repeat(1:5, inner = 5),
)
barplot!(ax4, tbl_gap.cat, tbl_gap.height,
    dodge = tbl_gap.grp,
    color = viscmap[tbl_gap.grp],
    direction = :x,
)

ax4.yticks = (eachindex(basins), basins)
# ax4.xticklabelrotation = pi/8
ax4.xlabel = "R-tipping gap (K)"
ax4.yaxisposition = :right
yl = (-0.6, 0)
# ylims!(ax4, yl)
text!(ax4, -0.5, -5, text = "c", font = :bold, color = :gray10, fontsize = 30)

labels = vcat(visc_labels...)
elements = [PolyElement(polycolor = viscmap[i]) for i in 1:5]
Legend(fig[0, 3], elements, labels, nbanks = 5, valign = :bottom)

cs = 470
rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -140)
rowgap!(fig.layout, 2, -170)
rowsize!(fig.layout, 0, 10)
rowsize!(fig.layout, 2, 10)
colsize!(fig.layout, 1, cs)
colsize!(fig.layout, 2, cs)
colsize!(fig.layout, 3, cs)

for ax in (ax1, ax4)
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xminorgridvisible = false
    ax.yminorgridvisible = false
end
xlims!(ax4, -0.6, 0)
xlims!(ax3, -0.6, 0)
ylims!(ax4, 0.5, 5.5)
ax4.xticks = -0.5:0.1:0
ax4.yreversed = true

fig

save(plotsdir("v2/rtip/rtip.png"), fig)
save(plotsdir("v2/rtip/rtip.pdf"), fig)

barplot!(ax4, tbl_bgap.cat, tbl_bgap.height,
    dodge = tbl_bgap.grp,
    color = viscmap[tbl_bgap.grp],
    direction = :x,
    alpha = 0.5,
)
fig

save(plotsdir("v2/rtip/rtip-btip.png"), fig)
save(plotsdir("v2/rtip/rtip-btip.pdf"), fig)