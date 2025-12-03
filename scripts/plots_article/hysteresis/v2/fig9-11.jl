include("../../../intro.jl")

#####################################################
# Common variables accross hysteresis figures
#####################################################

T = Float32
polar_amplification = 1.8
f_to = 0.25
f2015 = 1.2

s = 400                     # stride
lw1, lw2, lw3 = 5, 5, 7     # line widths
ms1, ms2 = 15, 20           # marker sizes

ssp126_2100 = 2.0
ssp245_2100 = 3.0
ssp370_2100 = 4.4
ssp585_2100 = 5.4
ssp2100 = [ssp126_2100, ssp245_2100, ssp370_2100, ssp585_2100]
ssp_labels = ["SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"]
ssp_line_opts = (linewidth = lw1, linestyle = :dash, alpha = 0.7)
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
    axx.xticks = -10:1:10
    axx.xminorticks = -10:0.2:20
    axx.xminorgridvisible = true
    axx.xlabel = L"GMT anomaly $f$ (K)"
    axx.yticks = 0:10:60
    axx.yminorticks = IntervalsBetween(10)
    axx.yminorgridvisible = true
    axx.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
    ylims!(axx, 0, 60)
    axislegend(axx, position = :rt)

    atm_axx.xticks = -10:2:30
    atm_axx.xgridvisible = false
    atm_axx.xaxisposition = :top
    atm_axx.xlabel = L"Regional atmospheric temperature anomaly $f_a$ (K)"
    hideydecorations!(atm_axx)
end

#####################################################
# Fig 9 - Comparison among processes
######################################################

dir = datadir("output/ais/v2/hyster")
xps = [
    "$dir/retreat/aqef/minvisc/how",
    "$dir/regrowth/aqef/how",
    "$dir/retreat/aqef/minvisc/dpr",
    "$dir/regrowth/aqef/dpr",
    "$dir/retreat/aqef/minvisc/atm",
    "$dir/regrowth/aqef/atm",
    "$dir/retreat/aqef/minvisc/refnomslow",
    "$dir/regrowth/aqef/refnomslow",
]
aqef = AQEFResults(T, xps)

lws = [lw2, lw2, lw2, lw2, lw2, lw2, lw3, lw3]
xp_labels = [
    "HOW",
    nothing,
    "DPR",
    nothing,
    "ATM",
    nothing,
    "REF",
    nothing,
]
cycling_colors = [
    xpcolors["HOW"],
    xpcolors["HOW"],
    xpcolors["DPR"],
    xpcolors["DPR"],
    xpcolors["ATM"],
    xpcolors["ATM"],
    xpcolors["REF"],
    xpcolors["REF"],
]

eqldir1 = "$dir/retreat/equil/refnomslow"
eql1 = EquilResults(T, eqldir1)
# eqldir2 = datadir("output/ais/hyster/16km/regrowth/equil")
# eql2 = EquilResults(T, eqldir2)

set_theme!(theme_latexfonts())
fig9 = Figure(size=(800, 800), fontsize = 22)
atm_ax1 = Axis(fig9[1, 1], aspect = AxisAspect(1))
ax1 = Axis(fig9[1, 1], aspect = AxisAspect(1))
for (i, ssp) in enumerate(ssp2100)
    vlines!(ax1, ssp, label = ssp_labels[i]*", 2100", color = xpcolors[ssp_labels[i]];
        ssp_line_opts...)
    # vlines!(ax1, ssp, color = xpcolors[ssp_labels[i]]; ssp_line_opts...)
end
for k in 1:aqef.n_xps
    lines!(ax1, aqef.f[k][1:s:end] ./ polar_amplification .+ f2015,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]))
end
x1, x2 = 0, 11
format_axs!(ax1, atm_ax1, x1, x2)
fig9
save(plotsdir("v2/hysteresis/fig9.png"), fig9)
save(plotsdir("v2/hysteresis/fig9.pdf"), fig9)

######################################################
# Fig 10 - Comparison to previous studies
######################################################

fig10 = Figure(size=(800, 800), fontsize = 22)
ax2 = Axis(fig10[1, 1], aspect = AxisAspect(1))
atm_ax2 = Axis(fig10[1, 1], aspect = AxisAspect(1))

h94dir = datadir("processed/huybrechts1994")
h94_retreat, _ = readdlm("$h94dir/h94-retreat.csv", ',', header = true)
h94_regrowth, _ = readdlm("$h94dir/h94-regrowth.csv", ',', header = true)
permindices = sortperm(h94_retreat[:, 1])
h94_retreat = h94_retreat[permindices, :]
permindices = sortperm(h94_regrowth[:, 1])
h94_regrowth = h94_regrowth[permindices, :]

g20dir = datadir("processed/garbe2020")
g20_retreat_ramp, _ = readdlm("$g20dir/retreat-ramp.csv", ',', header = true)
g20_retreat_equil, _ = readdlm("$g20dir/retreat-equil.csv", ',', header = true)
g20_regrowth_ramp, _ = readdlm("$g20dir/regrowth-ramp.csv", ',', header = true)
g20_regrowth_equil, _ = readdlm("$g20dir/regrowth-equil.csv", ',', header = true)

g20_retreat_equil = vcat([0, 55]', g20_retreat_equil)
g20_regrowth_equil = vcat([-3 55.6; -2 55.6; -1 55.3; 0 55], g20_regrowth_equil)

g20 = [g20_retreat_ramp, g20_retreat_equil, g20_regrowth_ramp, g20_regrowth_equil]
g20_labels = ["G20R", "G20E", nothing, nothing]
g20_colors = [xpcolors["G20R"], xpcolors["G20E"], xpcolors["G20R"], xpcolors["G20E"]]
g20_plotstyle = [:lines, :scatterlines, :lines, :scatterlines]
g20_markers = [nothing, :dtriangle, nothing, :utriangle]

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

m3_msle = 1e9 * 1e6 / (3.625 * 1e14)
# h94[1, 2] * m3_msle * gr_correction = 58
gr_correction = 58 / (h94_retreat[1, 2] * m3_msle)

scatterlines!(ax2, h94_retreat[:, 1] ./ polar_amplification,
    h94_retreat[:, 2] .* m3_msle .* gr_correction, marker = :dtriangle,
    linewidth = lw1, label = "H94", color = xpcolors["H94"], markersize = ms2)
scatterlines!(ax2, h94_regrowth[:, 1] ./ polar_amplification,
    h94_regrowth[:, 2] .* m3_msle .* gr_correction, marker = :utriangle,
    linewidth = lw1, color = xpcolors["H94"], markersize = ms2)



for k in [1, 2, 7, 8]
    lines!(ax2, aqef.f[k][1:s:end] ./ polar_amplification .+ f2015,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]))
end

# for axx in [ax2]
#     scatter!(axx, eql1.f ./ polar_amplification .+ f2015, eql1.V_sle;
#         color = :black, label = "EQL", markersize = ms1)
#     scatter!(axx, eql2.f ./ polar_amplification .+ f2015, eql2.V_sle;
#         color = :black, markersize = ms1)
# end

format_axs!(ax2, atm_ax2, x1, 12)
fig10
save(plotsdir("v2/hysteresis/fig10.png"), fig10)
save(plotsdir("v2/hysteresis/fig10.pdf"), fig10)

###############################################
# Fig 11 - Intermediate regrowth
###############################################

xps = [
    datadir("output/ais/v2/hyster/regrowth/aqef/intermediate/2"),
    datadir("output/ais/v2/hyster/retreat/aqef/minvisc/refnomslow"),
]
xp_labels = [
    "2",
    "REF",
]
aqef_intermediate = AQEFResults(T, xps)

lws = [lw2, lw3]
cycling_colors = [
    2,
    xpcolors["REF"],
]

fig11 = Figure(size=(800, 800), fontsize = 22)
atm_ax3 = Axis(fig11[1, 1], aspect = AxisAspect(1))
ax3 = Axis(fig11[1, 1], aspect = AxisAspect(1))

sref = 500
f_retreat = aqef_intermediate.f[end][1:sref:end] ./ polar_amplification .+ f2015
V_retreat = aqef_intermediate.V_sle[end][1:sref:end]
lines!(
    ax3,
    f_retreat,
    V_retreat,
    linewidth = lws[end],
    label = xp_labels[end],
    color = f_retreat,
    colormap = :jet,
    colorrange = (1, 11),
)
for i in 1:aqef_intermediate.n_xps-1
    f = aqef_intermediate.f[i][1:s:end] ./ polar_amplification .+ f2015
    V = aqef_intermediate.V_sle[i][1:s:end]
    c = fill(maximum(f), length(f))
    lines!(
        ax3,
        f,
        V,
        linewidth = lws[i],
        color = c,
        colormap = :jet,
        colorrange = (1, 11),
    )
end

format_axs!(ax3, atm_ax3, x1, x2)
fig11
save(plotsdir("v2/hysteresis/fig11.png"), fig11)
save(plotsdir("v2/hysteresis/fig11.pdf"), fig11)