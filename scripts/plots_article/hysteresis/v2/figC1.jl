include("../../../intro.jl")

# Metaparameters
T = Float32
polar_amplification = 1.8
f_to = 0.25
f2020 = 1.2
lw1, lw2 = 2, 5

# Load retreat
xps1 = [datadir("output/ais/v2/hyster/retreat/aqef/minvisc/refnomslow")]
aqef1 = AQEFResults(T, xps1)
s1 = Int(mean(diff(aqef1.t_2D[1])) / mean(diff(aqef1.t_1D[1])))
f1 = T.(aqef1.f[1][1:s1:end] ./ polar_amplification .+ f2020)

# Load regrowth
xps2 = [datadir("output/ais/v2/hyster/regrowth/aqef/refnomslow")]
aqef2 = AQEFResults(T, xps2)
s2 = Int(mean(diff(aqef2.t_2D[1])) / mean(diff(aqef2.t_1D[1])))
f2 = T.(aqef2.f[1][1:s2:end] ./ polar_amplification .+ f2020)

# Compute area-perimeter ratio
apr1 = get_area_perimeter_ratio(aqef1)
apr2 = get_area_perimeter_ratio(aqef2)

# Compute rate-velocity-slope product
grline1, dhvsp1 = rate_velocity_slope_product(xps1[1])
grline2, dhvsp2 = rate_velocity_slope_product(xps2[1])

# Post-process metrics. Smooth apr1 with 1:2:end because rate of forcing change is 2x faster than in regrowth
spline_level = 0.0001
apr1_smooth = fit(SmoothingSpline, Float64.(f1[1:2:end]), Float64.(apr1[1][1:2:end]), spline_level)
apr1_smooth = predict(apr1_smooth)
apr2_smooth = fit(SmoothingSpline, Float64.(f2), Float64.(apr2[1]), spline_level)
apr2_smooth = predict(apr2_smooth)

dapr1 = fill(NaN, length(apr1_smooth))
dapr2 = fill(NaN, length(apr2_smooth))
dapr1[2:end-1] .= (apr1_smooth[3:end] .- apr1_smooth[1:end-2]) / (2*dt1)
dapr2[2:end-1] .= (apr2_smooth[3:end] .- apr2_smooth[1:end-2]) / (dt2)

# fighh = Figure()
# axx1 = Axis(fighh[1, 1])
# axx2 = Axis(fighh[1, 2])
# axx3 = Axis(fighh[2, 1])
# axx4 = Axis(fighh[2, 2])
# lines!(axx1, f1, apr1[1])
# lines!(axx2, f2, apr2[1])
# lines!(axx1, f1, apr1_smooth)
# lines!(axx2, f2, apr2_smooth)
# lines!(axx3, f1[1:2:end], dapr1)
# lines!(axx4, f2, dapr2)
# ylims!(axx3, -2f-3, 2f-3)
# ylims!(axx4, -2f-3, 2f-3)
# fighh

dhvsp1_ts = MISIMetrics(aqef1.t_2D[1])
update!(dhvsp1_ts, dhvsp1)
dhvsp2_ts = MISIMetrics(aqef2.t_2D[1])
update!(dhvsp2_ts, dhvsp2)

cat_dapr = vcat(
    dapr1[not.(isnan.(dapr1)) .&& dapr1 .< 0],
    dapr2[not.(isnan.(dapr2)) .&& dapr2 .> 0])
cat_dhvsp = vcat(
    dhvsp1_ts.mini[not.(isnan.(dhvsp1_ts.mini))],
    dhvsp2_ts.maxi[not.(isnan.(dhvsp2_ts.maxi))])
pp_slp = 3 * std(cat_dhvsp)
pp_apr = 3 * std(cat_dapr)

####################################################################
# Plotting
#####################################################################

cmap = cgrad(:tab10, categorical = true)
slp_color = cmap[2]
apr_color = cmap[3]

set_theme!(theme_latexfonts())
figC1 = Figure(size=(1000, 800), fontsize = 16)
ax = Axis(figC1[1:2, 1])
ax_dslp = Axis(figC1[3, 1])
ax_dapr = Axis(figC1[4, 1])
ax2 = Axis(figC1[1:2, 2])
ax2_dslp = Axis(figC1[3, 2])
ax2_dapr = Axis(figC1[4, 2])

for axx in [ax, ax2, ax_dslp, ax2_dapr, ax_dapr, ax2_dslp]
    xlims!(axx, 0, 11)
    axx.xticks = -2:1:10
    axx.xminorticks = -2:0.2:12
    axx.xminorgridvisible = true
end

ax.title = "Retreat branch"
ax2.title = "Regrowth branch"
ax.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
for axx in [ax, ax2]
    axx.yticks = 10:10:60
    axx.yminorticks = IntervalsBetween(10)
    axx.yminorgridvisible = true
    axx.xaxisposition = :top
    ylims!(axx, 0, 60)
end
[axx.yaxisposition = :right for axx in [ax2, ax2_dapr, ax2_dslp]]

ax_dslp.ylabel = L"MISI metric $\varphi$ ($\mathrm{m^2 \, yr^{-2}}$)"
ax_dapr.ylabel = L"Perimeter metric $\psi$ ($\mathrm{yr^{-1}}$)"

for axx in [ax_dapr, ax2_dapr]
    axx.xlabel = L"GMT anomaly $f$ (K)"
end

for axx in [ax_dapr, ax2_dapr]
    xlims!(axx, 0, 11)
    ylims!(axx, -1.5f-3, 1.5f-3)
end
for axx in [ax_dslp, ax2_dslp]
    xlims!(axx, 0, 11)
    ylims!(axx, -15, 6)
    axx.xticklabelsvisible = false
    axx.xticksvisible = false
end

lgray = :gray75
dgray = :gray50
large_misi = [1.95, 5.87, 6.95, 8.65]
small_misi = [4.2, 6.21, 7.2, 7.6, 7.9, 9.2]
perimeter = [9.7, 10.1]
shadings = vcat(large_misi, small_misi, perimeter)


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
    (1.174, "Wilkes", 30, 7.5),
]
f_bif = [bifs[i][1] for i in eachindex(bifs)]
vlines!(ax, shadings, linewidth = 5, color = :gray60, alpha = 0.5)
vlines!(ax_dapr, shadings, linewidth = 5, color = :gray60, alpha = 0.5)
vlines!(ax_dslp, shadings, linewidth = 5, color = :gray60, alpha = 0.5)
vlines!(ax2, f_bif, color = :gray60, alpha = 0.5, linewidth = 5)
vlines!(ax2_dapr, f_bif, linewidth = 5, color = :gray60, alpha = 0.5)
vlines!(ax2_dslp, f_bif, linewidth = 5, color = :gray60, alpha = 0.5)

lines!(ax, f1, aqef1.V_sle[1][1:s1:end], linewidth = lw2, color = xpcolors["REF"])
lines!(ax2, f2, aqef2.V_sle[1][1:s2:end], linewidth = lw2, color = xpcolors["REF"])
# vlines!(ax, 1f6, color = dgray, linewidth = 5, alpha = 0.6, label = "Bifurcation")
# vlines!(ax, 1f6, color = lgray, linewidth = 5, alpha = 0.6, label = "Small Bifurcation")
# axislegend(ax, position = :lb)

alpha_band = 0.4
band!(ax_dslp, [0, 50], [-pp_slp, -pp_slp], [pp_slp, pp_slp], alpha = alpha_band, color = slp_color)
band!(ax_dapr, [0, 50], [-pp_apr, -pp_apr], [pp_apr, pp_apr], alpha = alpha_band, color = apr_color)
lines!(ax_dslp, f1, dhvsp1_ts.mini, linewidth = lw1, color = slp_color)
lines!(ax_dapr, f1[1:2:end], dapr1, linewidth = lw1, color = apr_color)

band!(ax2_dslp, [0, 50], [-pp_slp, -pp_slp], [pp_slp, pp_slp], alpha = alpha_band, color = slp_color)
band!(ax2_dapr, [0, 50], [-pp_apr, -pp_apr], [pp_apr, pp_apr], alpha = alpha_band, color = apr_color)
lines!(ax2_dslp, f2, dhvsp2_ts.maxi, linewidth = lw1, color = slp_color)
lines!(ax2_dapr, f2, dapr2, linewidth = lw1, color = apr_color)

[hlines!(axx, 0, color = :black, linestyle = :dash, linewidth = 1) for
    axx in [ax_dapr, ax_dslp, ax2_dapr, ax2_dslp]]

text!(ax, 0.2, 2, text = "(a)", font = :bold)
text!(ax2, 0.2, 2, text = "(b)", font = :bold)
text!(ax_dslp, 0.2, -14, text = "(c)", font = :bold)
text!(ax2_dslp, 0.2, -14, text = "(d)", font = :bold)
text!(ax_dapr, 0.2, -1.3f-3, text = "(e)", font = :bold)
text!(ax2_dapr, 0.2, -1.3f-3, text = "(f)", font = :bold)

rowgap!(figC1.layout, 5)
colgap!(figC1.layout, 5)
figC1

save(plotsdir("v2/hysteresis/figC1.png"), figC1)
save(plotsdir("v2/hysteresis/figC1.pdf"), figC1)