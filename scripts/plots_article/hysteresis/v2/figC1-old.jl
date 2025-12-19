include("../../../intro.jl")

T = Float32
polar_amplification = 1.8
f_to = 0.25
f2020 = 1.2
lw1, lw2 = 2, 5

xps = [datadir("output/ais/v2/hyster/retreat/aqef/minvisc/refnomslow")]
xp_labels = ["REF"]
aqef = AQEFResults(T, xps)

xps2 = [
    datadir("output/ais/v2/hyster/regrowth/aqef/refnomslow"),
    # datadir("output/ais/v2/hyster/regrowth/aqef/minvisc/refnomslow-restarted"),
]
xp_labels = [
    # nothing,
    "REF",
]
aqef2 = AQEFResults(T, xps2)
s = Int(mean(diff(aqef.t_2D[1])) / mean(diff(aqef.t_1D[1])))
# s2 = Int(mean(diff(aqef2.t_2D[1])) / mean(diff(aqef2.t_1D[1])))
ff = T.(aqef.f[1] ./ polar_amplification .+ f2020)
ff21 = T.(aqef2.f[1] ./ polar_amplification .+ f2020)
# ff22 = T.(aqef2.f[2] ./ polar_amplification .+ f2020)

apr = get_area_perimeter_ratio(aqef)
grline, flslp, zb_dx, zb_dy = get_slope_flow_product(xps[1], dslp = 1)
grline, asp = acceleration_slope_product(xps[1], dslp = 1)
grline, vasp = velocity_acceleration_slope_product(xps[1], dslp = 1)
k1 = 2
heatmap(view(zb_dx, :, :, k1))
heatmap(view(zb_dy, :, :, k1))
heatmap(view(flslp, :, :, k1))
heatmap(view(asp, :, :, k1))
heatmap(view(vasp, :, :, k1))
extrema(view(zb_dx, :, :, k1))
extrema(view(zb_dy, :, :, k1))
extrema(view(flslp, :, :, k1))
extrema(view(asp, :, :, k1))
extrema(view(vasp, :, :, k1)) 

k2 = 55
heatmap(view(zb_dx, :, :, k2))
heatmap(view(zb_dy, :, :, k2))
heatmap(view(flslp, :, :, k2))
heatmap(view(asp, :, :, k2))
heatmap(view(vasp, :, :, k2))
extrema(view(zb_dx, :, :, k2))
extrema(view(zb_dy, :, :, k2))
extrema(view(flslp, :, :, k2))
extrema(view(asp, :, :, k2))
extrema(view(vasp, :, :, k2))

mean_flslp_positive = Vector{Float32}(undef, length(aqef.t_2D[1]))
pnorm_flslp = Vector{Float32}(undef, length(aqef.t_2D[1]))
max_flslp = Vector{Float32}(undef, length(aqef.t_2D[1]))
percentile_flslp = Vector{Float32}(undef, length(aqef.t_2D[1]))
percentile_mean_flslp = Vector{Float32}(undef, length(aqef.t_2D[1]))
flslp_positive = max.(flslp, 0f0)
pnorm = 2

max_asp = Vector{Float32}(undef, length(aqef.t_2D[1]))
max_vasp = Vector{Float32}(undef, length(aqef.t_2D[1]))
min_vasp = Vector{Float32}(undef, length(aqef.t_2D[1]))
for k in axes(flslp, 3)
    mean_flslp_positive[k] = mean(view(flslp_positive, :, :, k))
    pnorm_flslp[k] = (mean(view(flslp_positive, :, :, k) .^ pnorm)) ^ (1 / pnorm)
    max_flslp[k] = maximum(view(flslp, :, :, k))
    percentile_flslp[k] = percentile(flslp[:, :, k][view(grline, :, :, k)], 99.5)
    percentile_mean_flslp[k] = mean(flslp[:, :, k][view(flslp, :, :, k) .> percentile_flslp[k]])

    max_asp[k] = maximum(view(asp, :, :, k))
    max_vasp[k] = maximum(view(vasp, :, :, k))
    min_vasp[k] = minimum(view(vasp, :, :, k))
end
figh = Figure()
ax1 = Axis(figh[1, 1])
ax2 = Axis(figh[2, 1])
ax3 = Axis(figh[3, 1])
lines!(ax1, aqef.f[1][1:s:end], max_flslp, label = "max")
lines!(ax1, aqef.f[1][1:s:end], percentile_flslp, label = "99.5 percentile")
lines!(ax1, aqef.f[1][1:s:end], percentile_mean_flslp, label = "mean above 99.5 percentile")
lines!(ax2, aqef.f[1][1:s:end], mean_flslp_positive, label = "mean positive")
lines!(ax3, aqef.f[1][1:s:end], pnorm_flslp, label = "pnorm")
axislegend(ax1)
axislegend(ax2)
axislegend(ax3)
figh
lines(max_asp)
lines(max_vasp)
lines(min_vasp)


apr2 = get_area_perimeter_ratio(aqef2)
grline21, flslp21 = get_slope_flow_product(xps2[1])
# grline22, flslp22 = get_slope_flow_product(xps2[2])
max_flslp1 = Vector{Float32}(undef, length(aqef2.t_2D[1]))
for k in axes(flslp21, 3)
    max_flslp1[k] = maximum(view(flslp21, :, :, k))
end

maxslp2 = Vector{Float32}(undef, length(aqef2.t_2D[2]))
for k in axes(flslp22, 3)
    maxslp2[k] = maximum(view(flslp22, :, :, k))
end

xp_idx = aqef2.n_xps
f_ref = aqef2.f[end] ./ polar_amplification .+ f2020
stiching_forcing = -0.4
stiching_idx = findfirst(f_ref .<= stiching_forcing)



dt1 = mean(diff(aqef.t_2D[1]))
dt2 = mean(diff(aqef2.t_2D[1]))
smoothing = "loess"

if smoothing == "loess"
    dloess = 0.02
    mslp = loess(ff[1:s:end], ff[1:s:end], max_flslp, span = dloess)
    mslp1 = loess(ff21[1:s2:end], ff21[1:s2:end], max_flslp1, span = dloess)
    mslp2 = loess(ff22[1:s2:end], ff22[1:s2:end], maxslp2, span = dloess)
    mapr = loess(ff[1:s:end], ff[1:s:end], apr[1], span = dloess)
    mapr1 = loess(ff21[1:s2:end], ff21[1:s2:end], apr2[1], span = dloess)
    mapr2 = loess(ff22[1:s2:end], ff22[1:s2:end], apr2[2], span = dloess)
elseif smoothing == "rollmean"
    nx = 2
    mslp = rollmean(ff[1:s:end], max_flslp, nx = nx * 2)
    mslp1 = rollmean(ff21[1:s2:end], max_flslp1, nx = nx)
    mslp2 = rollmean(ff22[1:s2:end], maxslp2, nx = nx)
    mapr = rollmean(ff[1:s:end], apr[1], nx = nx * 2)
    mapr1 = rollmean(ff21[1:s2:end], apr2[1], nx = nx)
    mapr2 = rollmean(ff22[1:s2:end], apr2[2], nx = nx)
end

detrend = "symmeandiff_smooth"     # ["diff", "symdiff", "symdiff_smooth", "diff_smooth", "detrend"]

if detrend == "diff_smooth"
    dmslp = diff(mslp)
    dmslp1 = diff(mslp1)
    dmslp2 = diff(mslp2)
    dmapr = diff(mapr)
    dmapr1 = diff(mapr1)
    dmapr2 = diff(mapr2)
    k1, k2 = 2, s2
elseif detrend == "symdiff_smooth"
    n = 1
    dmslp = sym_diff(mslp, n, dt1)
    dmslp1 = sym_diff(mslp1, n, dt2)
    dmslp2 = sym_diff(mslp2, n, dt2)
    dmapr = sym_diff(mapr, n, dt1)
    dmapr1 = sym_diff(mapr1, n, dt2)
    dmapr2 = sym_diff(mapr2, n, dt2)
    k1, k2 = 1, 0
elseif detrend == "symmeandiff_smooth"
    n = 1
    dmslp = sym_mean_diff(mslp, 2*n, dt1)
    dmslp1 = sym_mean_diff(mslp1, n, dt2)
    dmslp2 = sym_mean_diff(mslp2, n, dt2)
    dmapr = sym_mean_diff(mapr, 2*n, dt1)
    dmapr1 = sym_mean_diff(mapr1, n, dt2)
    dmapr2 = sym_mean_diff(mapr2, n, dt2)
    k1, k2 = 1, 0
elseif detrend == "detrend"
    dmslp = max_flslp .- mslp
    dmslp1 = max_flslp1 .- mslp1
    dmslp2 = maxslp2 .- mslp2
    dmapr = apr[1] .- mapr
    dmapr1 = apr2[1] .- mapr1
    dmapr2 = apr2[2] .- mapr2
    k1, k2 = 1, 0
end

cat_dmslp = vcat(
    dmslp[not.(isnan.(dmslp))],
    dmslp1[not.(isnan.(dmslp1))],
    dmslp2[not.(isnan.(dmslp2))])
cat_dmapr = vcat(
    dmapr[not.(isnan.(dmapr))],
    dmapr1[not.(isnan.(dmapr1))],
    dmapr2[not.(isnan.(dmapr2))])
pp_slp = percentile(abs.(cat_dmslp), 98)
pp_apr = percentile(abs.(cat_dmapr), 98)
pp_slp = 3 * std(cat_dmslp)
pp_apr = 3 * std(cat_dmapr)

####################################################################
# Plotting
#####################################################################

set_theme!(theme_latexfonts())
color_scheme = cgrad(:seaborn_colorblind, 10, categorical = true)
apr_color = color_scheme[2]
slp_color = color_scheme[3]
alpha_band = 0.4
lw1, lw2 = 2, 4

figA4 = Figure(size=(1000, 800), fontsize = 20)
ax = Axis(figA4[1:2, 1])
ax_slp = Axis(figA4[3, 1])
ax_apr = Axis(figA4[3, 1])
ax_dslp = Axis(figA4[4, 1])
ax_dapr = Axis(figA4[4, 1])

ax2 = Axis(figA4[1:2, 2])
ax2_slp = Axis(figA4[3, 2])
ax2_apr = Axis(figA4[3, 2])
ax2_dapr = Axis(figA4[4, 2])
ax2_dslp = Axis(figA4[4, 2])

for axx in [ax, ax_slp, ax_apr, ax_dslp, ax_dapr]
    xlims!(axx, 0, 11)
    ax.xticks = -2:2:12
    ax.xminorticks = -2:0.2:12
end
for axx in [ax2, ax2_slp, ax2_apr, ax2_dapr, ax2_dslp]
    xlims!(axx, 0, 11)
    axx.xticks = -2:2:12
    ax.xminorticks = -2:0.2:12
end
[hidedecorations!(axx) for axx in [ax_apr, ax_dapr, ax2_slp, ax2_dslp]]
[ylims!(axx, (-0.5, 1.5) ./ 1f3) for axx in [ax2_dapr, ax_dapr]]
[ylims!(axx, (-3, 1.5) ./ 1f3) for axx in [ax2_dslp, ax_dslp]]
[ylims!(axx, (-1, 20)) for axx in [ax_slp, ax2_slp]]
[ylims!(axx, (-1, 45)) for axx in [ax_apr, ax2_apr]]
ax_dslp.yticks = -2f-3:2f-3:2f-3
ax2_dapr.yticks = -5f-4:5f-4:1f-3

ax.title = "Retreat branch"
ax2.title = "Regrowth branch"
for axx in [ax, ax2]
    axx.yticks = 0:10:60
    axx.yminorticks = IntervalsBetween(10)
    axx.xminorgridvisible = true
    axx.yminorgridvisible = true
    axx.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
    axx.xaxisposition = :top
    ylims!(axx, 0, 60)
end
[axx.yaxisposition = :right for axx in [ax2, ax2_slp, ax2_apr, ax2_dapr, ax2_dslp]]

ax_slp.ylabel = L"MISI metric $\varphi$ ($\mathrm{m \, yr^{-1}}$)"
ax_dslp.ylabel = L"$\dot{\varphi}$ ($\mathrm{m \, yr^{-2}}$)"
for axx in [ax_slp, ax_dslp]
    axx.yticklabelcolor = slp_color
    axx.ytickcolor = slp_color
    axx.ylabelcolor = slp_color
    axx.ygridvisible = false
end

ax2_apr.ylabel = L"Perimeter metric $\psi$ (1)"
ax2_dapr.ylabel = L"$\dot{\psi}$ ($\mathrm{yr^{-1}}$)"
for axx in [ax2_apr, ax2_dapr]
    axx.yticklabelcolor = apr_color
    axx.ytickcolor = apr_color
    axx.ylabelcolor = apr_color
    axx.ygridvisible = false
end

for axx in [ax_apr, ax_slp, ax2_apr, ax2_slp]
    axx.xticklabelsvisible = false
    axx.xticksvisible = false
end

for axx in [ax_dslp, ax2_dapr]
    axx.xlabel = L"GMT anomaly $f$ (K)"
end

lgray = :gray75
dgray = :gray50
# shade = [(1.2, 1.3), (4.7, 4.8), (6.1, 6.2), (7.0, 7.1), (7.8, 7.9)]
# lightshade = [(4.6, 4.7), (8.5, 8.6), (9.9, 10), (10.2, 10.3)]
# shade_transitions!([ax_slp, ax_dslp, ax], shade, f2020, 0.2, dgray)
# shade_transitions!([ax_slp, ax_dslp, ax], lightshade, 0.0, 0.2, lgray)

# shade = [(1.3, 1.4), (2.6, 2.7), (2.9, 3.0), (3.25, 3.35),
#     (3.6, 3.7), (4.1, 4.2), (4.8, 4.9), (8.05, 8.15), (8.65, 8.75)]
# lightshade = [(0, 0.1), (0.35, 0.45), (6.6, 6.7), (9.1, 9.2)]
# shade_transitions!([ax2_slp, ax2_dslp, ax2], shade, 0.0, 0.2, dgray)
# shade_transitions!([ax2_slp, ax2_dslp, ax2], lightshade, 0.0, 0.2, lgray)

lines!(ax2, aqef2.f[1][1:s2:end] ./ polar_amplification .+ f2020,
    aqef2.V_sle[1][1:s2:end], linewidth = lw2, label = xp_labels[1],
    color = xpcolors["REF"])
vlines!(ax, 1f6, color = dgray, linewidth = 5, alpha = 0.6, label = "Bifurcation")
vlines!(ax, 1f6, color = lgray, linewidth = 5, alpha = 0.6, label = "Small Bifurcation")
axislegend(ax, position = :lb)

lines!(ax, ff[1:s:end], aqef.V_sle[1][1:s:end], linewidth = lw2,
    label = xp_labels[1], color = xpcolors["REF"])
lines!(ax_slp, ff[1:s:end], max_flslp, linewidth = lw2, color = slp_color, alpha = 0.4)
lines!(ax_slp, ff[1:s:end], mslp, linewidth = lw1, color = slp_color)
lines!(ax_apr, ff[1:s:end][2:end], apr[1][2:end], linewidth = lw2, color = apr_color, alpha = 0.4)
lines!(ax_apr, ff[1:s:end], mapr, linewidth = lw1, color = apr_color)
band!(ax_dslp, [0, 50], [-pp_slp, -pp_slp], [pp_slp, pp_slp], alpha = alpha_band, color = slp_color)
band!(ax_dapr, [0, 50], [-pp_apr, -pp_apr], [pp_apr, pp_apr], alpha = alpha_band, color = apr_color)
lines!(ax_dapr, ff[k1:s:end], dmapr, linewidth = lw1, color = apr_color)
lines!(ax_dslp, ff[k1:s:end], dmslp, linewidth = lw1, color = slp_color)

lines!(ax2_slp, ff21[1:s2:end], max_flslp1, linewidth = lw2, color = slp_color, alpha = 0.4)
lines!(ax2_slp, ff22[1:s2:end], maxslp2, linewidth = lw2, color = slp_color, alpha = 0.4)
lines!(ax2_slp, ff21[1:s2:end], mslp1, linewidth = lw1, color = slp_color)
lines!(ax2_slp, ff22[1:s2:end], mslp2, linewidth = lw1, color = slp_color)
lines!(ax2_apr, ff21[1:s2:end], apr2[1], linewidth = lw2, color = apr_color, alpha = 0.4)
lines!(ax2_apr, ff22[1:s2:end], apr2[2], linewidth = lw2, color = apr_color, alpha = 0.4)
lines!(ax2_apr, ff22[1:s2:end], mapr2, linewidth = lw1, color = apr_color)
lines!(ax2_apr, ff21[1:s2:end], mapr1, linewidth = lw1, color = apr_color)

band!(ax2_dslp, [-50, 50], [-pp_slp, -pp_slp], [pp_slp, pp_slp],
    alpha = alpha_band, color = slp_color)
band!(ax2_dapr, [-50, 50], [-pp_apr, -pp_apr], [pp_apr, pp_apr],
    alpha = alpha_band, color = apr_color)
lines!(ax2_dslp, ff21[k1:s2:end], dmslp1, linewidth = lw1, color = slp_color)
lines!(ax2_dslp, ff22[k1:s2:end-k2], dmslp2, linewidth = lw1, color = slp_color)
lines!(ax2_dapr, ff21[k1:s2:end], dmapr1, linewidth = lw1, color = apr_color)
lines!(ax2_dapr, ff22[k1:s2:end-k2], dmapr2, linewidth = lw1, color = apr_color)
[hlines!(axx, 0, color = :black, linestyle = :dash, linewidth = 1) for
    axx in [ax_dapr, ax_dslp, ax2_dapr, ax2_dslp]]

figA4

save(plotsdir("16km/hysteresis/figA4.png"), figA4)
save(plotsdir("16km/hysteresis/figA4.pdf"), figA4)