include("../intro.jl")

T = Float32
xps = [
    datadir("output/ais/v2/hyster/retreat/aqef/atm"),
    datadir("output/ais/v2/hyster/retreat/aqef/ocn"),
    datadir("output/ais/v2/hyster/retreat/aqef/how"),
    datadir("output/ais/v2/hyster/retreat/aqef/dpr"),
    datadir("output/ais/v2/hyster/retreat/aqef/refm2"),
    datadir("output/ais/v2/hyster/retreat/aqef/refp2"),
    datadir("output/ais/v2/hyster/retreat/aqef/refp2s"),
    datadir("output/ais/v2/hyster/retreat/aqef/refnom"),
]

xp_labels = [
    "ATM",
    "OCN",
    "HOW",
    "DPR",
    L"$-2 \, \sigma$",
    L"$+2 \, \sigma$",
    L"$+2 \, \sigma$ slow",
    "REF",
]
aqef = AQEFResults(T, xps)

# eqldir1 = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
# eql1 = EquilResults(T, eqldir1)
# eqldir2 = datadir("output/ais/hyster/16km/regrowth/equil")
# eql2 = EquilResults(T, eqldir2)

lw1, lw2, lw3 = 5, 5, 7
ms1, ms2 = 15, 20
lws = [lw2, lw2, lw2, lw2, lw2, lw2, lw2, lw3]
cycling_colors = [
    xpcolors["ATM"],
    xpcolors["OCN"],
    xpcolors["HOW"],
    xpcolors["DPR"],
    xpcolors["REFm2"],
    xpcolors["REFp2"],
    :red,
    xpcolors["REF"],
]

polar_amplification = 1.8
f_to = 0.25
f_ref = aqef.f[end] ./ polar_amplification

ssp126_2100 = 2.0
ssp245_2100 = 3.0
ssp370_2100 = 4.4
ssp585_2100 = 5.4
ssp2100 = [ssp126_2100, ssp245_2100, ssp370_2100, ssp585_2100]
ssp_labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]

set_theme!(theme_latexfonts())
fig = Figure(size=(1400, 800), fontsize = 26)
# atm_ax1 = Axis(fig[1, 1], aspect = AxisAspect(1))
# hideydecorations!(atm_ax1)
ax1 = Axis(fig[1, 1], aspect = AxisAspect(1))
ax1.title = "(a) Comparison among processes"

# atm_ax2 = Axis(fig[1, 2], aspect = AxisAspect(1))
# hideydecorations!(atm_ax2)
# atm_ax2.xaxisposition = :top
# ax2 = Axis(fig[1, 2], aspect = AxisAspect(1))
# ax2.title = "(b) Comparison to previous studies"
# ax2.yaxisposition = :right

# ssp_line_opts = (linewidth = lw1, linestyle = :dash)
# for (i, ssp) in enumerate(ssp2100)
#     vlines!(ax1, ssp, label = ssp_labels[i], color = xpcolors[ssp_labels[i]];
#         ssp_line_opts...)
#     vlines!(ax1, ssp, color = xpcolors[ssp_labels[i]]; ssp_line_opts...)
# end

# h94dir = datadir("processed/huybrechts1994")
# h94_retreat, _ = readdlm("$h94dir/h94-retreat.csv", ',', header = true)
# h94_regrowth, _ = readdlm("$h94dir/h94-regrowth.csv", ',', header = true)
# permindices = sortperm(h94_retreat[:, 1])
# h94_retreat = h94_retreat[permindices, :]
# permindices = sortperm(h94_regrowth[:, 1])
# h94_regrowth = h94_regrowth[permindices, :]

# g20dir = datadir("processed/garbe2020")
# g20_retreat_ramp, _ = readdlm("$g20dir/retreat-ramp.csv", ',', header = true)
# g20_retreat_equil, _ = readdlm("$g20dir/retreat-equil.csv", ',', header = true)
# g20_regrowth_ramp, _ = readdlm("$g20dir/regrowth-ramp.csv", ',', header = true)
# g20_regrowth_equil, _ = readdlm("$g20dir/regrowth-equil.csv", ',', header = true)

# g20_retreat_equil = vcat([0, 55]', g20_retreat_equil)
# g20_regrowth_equil = vcat([-3 55.6; -2 55.6; -1 55.3; 0 55], g20_regrowth_equil)

# g20 = [g20_retreat_ramp, g20_retreat_equil, g20_regrowth_ramp, g20_regrowth_equil]
# g20_labels = ["G20R", "G20E", nothing, nothing]
# g20_colors = [xpcolors["G20R"], xpcolors["G20E"], xpcolors["G20R"], xpcolors["G20E"]]
# g20_plotstyle = [:lines, :scatterlines, :lines, :scatterlines]
# g20_markers = [nothing, :dtriangle, nothing, :utriangle]

# for i in eachindex(g20)
#     if g20_plotstyle[i] == :lines
#         lines!(ax2, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2], linewidth = 3,
#             label = g20_labels[i], color = g20_colors[i])
#     elseif g20_plotstyle[i] == :scatterlines
#         scatterlines!(ax2, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2],
#             linewidth = 3, label = g20_labels[i], color = g20_colors[i],
#             marker = g20_markers[i], markersize = ms2)
#     end
# end

# m3_msle = 1e9 * 1e6 / (3.625 * 1e14)
# # h94[1, 2] * m3_msle * gr_correction = 58
# gr_correction = 58 / (h94_retreat[1, 2] * m3_msle)

# scatterlines!(ax2, h94_retreat[:, 1] ./ polar_amplification,
#     h94_retreat[:, 2] .* m3_msle .* gr_correction, marker = :dtriangle,
#     linewidth = lw1, label = "H94", color = xpcolors["H94"], markersize = ms2)
# scatterlines!(ax2, h94_regrowth[:, 1] ./ polar_amplification,
#     h94_regrowth[:, 2] .* m3_msle .* gr_correction, marker = :utriangle,
#     linewidth = lw1, color = xpcolors["H94"], markersize = ms2)

s = 50
f2015 = 1.2
uniqueidx(v) = unique(i -> v[i], eachindex(v))
for k in 1:aqef.n_xps
    # idx = (aqef.f[k] .< 25) .&& (0 .< aqef.V_sle[k] .< 80)
    idx = uniqueidx(aqef.t_1D[k])
    lines!(ax1, aqef.f[k][idx][1:s:end] ./ polar_amplification .+ f2015,
        aqef.V_sle[k][idx][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]))
end

# for k in [1, 2, 7, 8, 9]
#     if k < aqef.n_xps
#         lines!(ax2, aqef.f[k][1:s:end] ./ polar_amplification .+ f2015,
#             aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
#             color = lcolor(cycling_colors[k]))
#     else
#         lines!(ax2, aqef.f[k][1:s:stiching_idx] ./ polar_amplification .+ f2015,
#             aqef.V_sle[k][1:s:stiching_idx], linewidth = lws[k], label = xp_labels[k],
#             color = lcolor(cycling_colors[k]))
#     end
# end

# for axx in [ax2]
#     scatter!(axx, eql1.f ./ polar_amplification .+ f2015, eql1.V_sle;
#         color = :black, label = "EQL", markersize = ms1)
#     scatter!(axx, eql2.f ./ polar_amplification .+ f2015, eql2.V_sle;
#         color = :black, markersize = ms1)
# end

# x1, x2 = 0, 12
# for axx in [ax1, ax2]
#     xlims!(axx, x1, x2)
# end
# for axx in [atm_ax1, atm_ax2]
#     xlims!(axx, x1 * polar_amplification, x2 * polar_amplification)
# end

# for axx in [ax1, ax2]
#     axx.titlegap = 80
#     axx.xticks = -10:2:10
#     axx.xminorticks = -10:0.2:20
#     axx.xminorgridvisible = true
#     axx.xlabel = L"GMT anomaly $f$ (K)"
#     axx.yticks = 0:10:60
#     axx.yminorticks = IntervalsBetween(10)
#     axx.yminorgridvisible = true
#     axx.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
#     ylims!(axx, 0, 60)
#     axislegend(axx, position = :rt)
# end

# for axx in [atm_ax1, atm_ax2]
#     axx.xticks = -10:2:30
#     axx.xgridvisible = false
#     axx.xaxisposition = :top
#     axx.xlabel = L"Regional atm. temperature anomaly $f_a$ (K)"
# end

axislegend(ax1, position = :lb)
xlims!(ax1, 0, 12)
ylims!(ax1, 0, 60)
fig
save(plotsdir("v2/16km/hysteresis/Haf_compare-hyst.png"), fig)
# save(plotsdir("16km/hysteresis/compare-hyst.pdf"), fig)