include("../../../intro.jl")

#####################################################
# Common variables accross hysteresis figures
#####################################################

T = Float32
polar_amplification = 1.8
f_to = 0.25
f2020 = 1.2

s = 400                     # stride
lw1, lw2, lw3 = 5, 5, 7     # line widths
ms1, ms2 = 8, 20           # marker sizes

# From tab spm1 of ipcc
ssp126_2100 = 1.8
ssp245_2100 = 2.7
ssp370_2100 = 3.6
ssp585_2100 = 4.4

ssp126_2100_hi = 2.4
ssp245_2100_hi = 3.5
ssp370_2100_hi = 4.6
ssp585_2100_hi = 5.7

ssp2100 = [ssp126_2100, ssp245_2100, ssp370_2100, ssp585_2100]
ssp2100_hi = [ssp126_2100_hi, ssp245_2100_hi,
    ssp370_2100_hi, ssp585_2100_hi]
ssp_labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
ssp_short_labels = [
    L"\bar{f}_\text{SSP1-2.6}",
    L"\bar{f}_\text{SSP2-4.5}",
    L"\bar{f}_\text{SSP3-7.0}",
    L"\bar{f}_\text{SSP5-8.5}",
]
ssphi_short_labels = [
    L"\hat{f}_\text{SSP1-2.6}",
    L"\hat{f}_\text{SSP2-4.5}",
    L"\hat{f}_\text{SSP3-7.0}",
    L"\hat{f}_\text{SSP5-8.5}",
]
ssp_line_opts = (linewidth = lw1, alpha = 0.7)
ssphi_line_opts = (linestyle = :dash, linewidth = 3, alpha = 0.7)
ssp_colors = [
    :darkblue,
    :orange,
    :red,
    :darkred,
]

function format_axs!(axx, atm_axx, x1, x2)

    xlims!(axx, x1, x2)
    xlims!(atm_axx, x1 * polar_amplification, x2 * polar_amplification)

    axx.titlegap = 80
    axx.xticks = -10:1:12
    axx.xminorticks = -10:0.2:20
    axx.xminorgridvisible = true
    axx.xlabel = L"GMT anomaly $f$ (K)"
    axx.yticks = 0:10:60
    axx.yminorticks = IntervalsBetween(10)
    axx.yminorgridvisible = true
    axx.ylabel = L"Equilibrium AIS volume $V_\mathrm{af}$ (m SLE)"
    ylims!(axx, 0, 60)

    atm_axx.xticks = -10:2:30
    atm_axx.xgridvisible = false
    atm_axx.xaxisposition = :top
    atm_axx.xlabel = L"Regional atmospheric temperature anomaly $f_a$ (K)"
    hideydecorations!(atm_axx)
end

#####################################################
# Loading data
######################################################

dir = datadir("output/ais/v2/hyster")
xps = [
    "$dir/retreat/aqef/refnomslow",
    "$dir/regrowth/aqef/refnomslow",
    "$dir/regrowth/aqef/refnomslow-restarted",
]
aqef = AQEFResults(T, xps)

lws = [lw3, lw3, lw3]
xp_labels = [
    "REF",
    nothing,
    nothing,
]
cycling_colors = [
    xpcolors["REF"],
    xpcolors["REF"],
    xpcolors["REF"],
]

eqldir1 = "$dir/retreat/equil/refnomslow"
eql1 = EquilResults(T, eqldir1)
eqldir2 = datadir("output/ais/v2/hyster/regrowth/equil/refnomslow")
eql2 = EquilResults(T, eqldir2)

######################################################
# Fig 11 - Comparison to previous studies
######################################################

set_theme!(theme_latexfonts())
x1, x2 = 0, 11
fig11 = Figure(size=(800, 800), fontsize = 22)
ax2 = Axis(fig11[1, 1], aspect = AxisAspect(1))
atm_ax2 = Axis(fig11[1, 1], aspect = AxisAspect(1))

g20dir = datadir("processed/garbe2020")
g20_retreat_ramp, _ = readdlm("$g20dir/retreat-ramp.csv", ',', header = true)
g20_retreat_equil, _ = readdlm("$g20dir/retreat-equil.csv", ',', header = true)
g20_regrowth_ramp, _ = readdlm("$g20dir/regrowth-ramp.csv", ',', header = true)
g20_regrowth_equil, _ = readdlm("$g20dir/regrowth-equil.csv", ',', header = true)

w26_retreat_ramp, _ = readdlm(datadir("processed/winkelmann2026/retreat.csv"), ',', header = false)

g20_retreat_equil = vcat([0, 55]', g20_retreat_equil)
g20_regrowth_equil = vcat([-3 55.6; -2 55.6; -1 55.3; 0 55], g20_regrowth_equil)

# g20 = [g20_retreat_equil, g20_regrowth_equil]
# g20_labels = ["G20E", nothing]
# g20_colors = [xpcolors["G20E"], xpcolors["G20E"]]
# g20_plotstyle = [:scatterlines, :scatterlines]
# g20_markers = [:dtriangle, :utriangle]
g20 = [g20_retreat_ramp, g20_retreat_equil, g20_regrowth_ramp, g20_regrowth_equil]
g20_labels = ["G20R", "G20E", nothing, nothing]
g20_colors = [xpcolors["G20R"], xpcolors["G20E"], xpcolors["G20R"], xpcolors["G20E"]]
g20_plotstyle = [:lines, :scatterlines, :lines, :scatterlines]
g20_markers = [nothing, :dtriangle, nothing, :utriangle]

for (i, ssp) in enumerate(ssp2100)
    vlines!(ax2, ssp, label = ssp_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssp_line_opts...)
    vlines!(ax2, ssp2100_hi[i], label = ssphi_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssphi_line_opts...)
    # vlines!(ax1, ssp, color = xpcolors[ssp_labels[i]]; ssp_line_opts...)
end

for i in eachindex(g20)
    if g20_plotstyle[i] == :lines
        lines!(ax2, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2], linewidth = 3,
            label = g20_labels[i], color = g20_colors[i])
    elseif g20_plotstyle[i] == :scatterlines
        scatterlines!(ax2, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2],
            linewidth = 3, label = g20_labels[i], color = g20_colors[i],
            marker = g20_markers[i], markersize = ms2)
    end
end

# plot_winkelmann2026 = true
# if plot_winkelmann2026
#     lines!(ax2, w26_retreat_ramp[:, 1] ./ polar_amplification, w26_retreat_ramp[:, 2], linewidth = 3,
#         label = "W26R", color = xpcolors["W26"])
# end

m3_msle = 1e9 * 1e6 / (3.625 * 1e14)
# h94[1, 2] * m3_msle * gr_correction = 58
# gr_correction = 58 / (h94_retreat[1, 2] * m3_msle)

for k in 1:3
    lines!(ax2, aqef.f[k][1:s:end] ./ polar_amplification .+ f2020,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]))
end

idx1 = sortperm(eql1.f)
f1 = eql1.f[idx1]
V_sle1 = eql1.V_sle[idx1]
V_sle1 = filter_outliers(V_sle1, mode = :low)
scatter!(ax2, f1 ./ polar_amplification .+ f2020, V_sle1;
    color = xpcolors["EQL"], label = "EQL", markersize = ms1)
idx2 = sortperm(eql2.f)
f2 = eql2.f[idx2]
V_sle2 = eql2.V_sle[idx2]
V_sle2 = filter_outliers(V_sle2, mode = :low)
scatter!(ax2, f2 ./ polar_amplification .+ f2020, V_sle2;
    color = :cornflowerblue, markersize = ms1)

format_axs!(ax2, atm_ax2, x1, x2)
axislegend(ax2, position = :rt, nbanks = 2)

fig11
save(plotsdir("v2/hysteresis/fig11-reduced.png"), fig11)
save(plotsdir("v2/hysteresis/fig11-reduced.pdf"), fig11)