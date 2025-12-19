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
    axx.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
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
    "$dir/retreat/aqef/minvisc/how",
    "$dir/regrowth/aqef/how",
    "$dir/retreat/aqef/minvisc/dpr",
    "$dir/regrowth/aqef/dpr",
    "$dir/retreat/aqef/minvisc/atm",
    "$dir/regrowth/aqef/atm",
    "$dir/retreat/aqef/minvisc/refnomslow",
    "$dir/regrowth/aqef/refnomslow",
    "$dir/regrowth/aqef/minvisc/refnomslow-restarted",
]
aqef = AQEFResults(T, xps)

lws = [lw2, lw2, lw2, lw2, lw2, lw2, lw3, lw3, lw3]
xp_labels = [
    "HOW",
    nothing,
    "DPR",
    nothing,
    "ATM",
    nothing,
    "REF",
    nothing,
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
    xpcolors["REF"],
]

eqldir1 = "$dir/retreat/equil/refnomslow"
eql1 = EquilResults(T, eqldir1)
eqldir2 = datadir("output/ais/v2/hyster/regrowth/equil/refnomslow")
eql2 = EquilResults(T, eqldir2)

#####################################################
# Fig 10 - Comparison among processes
######################################################

set_theme!(theme_latexfonts())
fig10 = Figure(size=(800, 800), fontsize = 22)
atm_ax1 = Axis(fig10[1, 1], aspect = AxisAspect(1))
ax1 = Axis(fig10[1, 1], aspect = AxisAspect(1))
for (i, ssp) in enumerate(ssp2100)
    vlines!(ax1, ssp, label = ssp_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssp_line_opts...)
    vlines!(ax1, ssp2100_hi[i], label = ssphi_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssphi_line_opts...)
    # vlines!(ax1, ssp, color = xpcolors[ssp_labels[i]]; ssp_line_opts...)
end
for k in 3:aqef.n_xps
    lines!(ax1, aqef.f[k][1:s:end] ./ polar_amplification .+ f2020,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]))
end
x1, x2 = 0, 11
format_axs!(ax1, atm_ax1, x1, x2)
axislegend(ax1, position = :rt)

fig10
save(plotsdir("v2/hysteresis/fig10.png"), fig10)
save(plotsdir("v2/hysteresis/fig10.pdf"), fig10)

######################################################
# Fig 11 - Comparison to previous studies
######################################################

fig11 = Figure(size=(800, 800), fontsize = 22)
ax2 = Axis(fig11[1, 1], aspect = AxisAspect(1))
atm_ax2 = Axis(fig11[1, 1], aspect = AxisAspect(1))

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

m3_msle = 1e9 * 1e6 / (3.625 * 1e14)
# h94[1, 2] * m3_msle * gr_correction = 58
gr_correction = 58 / (h94_retreat[1, 2] * m3_msle)

scatterlines!(ax2, h94_retreat[:, 1] ./ polar_amplification,
    h94_retreat[:, 2] .* m3_msle .* gr_correction, marker = :dtriangle,
    linewidth = lw1, label = "H94", color = xpcolors["H94"], markersize = ms2)
scatterlines!(ax2, h94_regrowth[:, 1] ./ polar_amplification,
    h94_regrowth[:, 2] .* m3_msle .* gr_correction, marker = :utriangle,
    linewidth = lw1, color = xpcolors["H94"], markersize = ms2)


for k in [1, 2, 7, 8, 9]
    lines!(ax2, aqef.f[k][1:s:end] ./ polar_amplification .+ f2020,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]))
end

idx1 = sortperm(eql1.f)
f1 = eql1.f[idx1]
V_sle1 = eql1.V_sle[idx1]
V_sle1 = filter_outliers(V_sle1)
scatter!(ax2, f1 ./ polar_amplification .+ f2020, V_sle1;
    color = xpcolors["EQL"], label = "EQL", markersize = ms1)
idx2 = sortperm(eql2.f)
f2 = eql2.f[idx2]
V_sle2 = eql2.V_sle[idx2]
V_sle2 = filter_outliers(V_sle2, mode = :low)
scatter!(ax2, f2 ./ polar_amplification .+ f2020, V_sle2;
    color = :cornflowerblue, markersize = ms1)

format_axs!(ax2, atm_ax2, x1, 12)
axislegend(ax2, position = :rt, nbanks = 2)

fig11
save(plotsdir("v2/hysteresis/fig11.png"), fig11)
save(plotsdir("v2/hysteresis/fig11.pdf"), fig11)

###############################################
# Fig 11 - Intermediate regrowth
###############################################

xps = readdir(datadir("output/ais/v2/hyster/regrowth/aqef/minvisc/intermediate"),
    join=true)

function extract_xp_number(str)
    parts = split(str, '/')
    lastpart = parts[end]
    return parse(Int, last(split(lastpart, '_')))
end

ii = extract_xp_number.(xps)
iii = sortperm(ii)
xps = xps[iii]
push!(xps, datadir("output/ais/v2/hyster/retreat/aqef/minvisc/refnomslow"))
aqef_intermediate = AQEFResults(T, xps)
lws = fill(lw2, length(ii))
push!(lws, lw3)

fig12 = Figure(size=(800, 800), fontsize = 22)
atm_ax3 = Axis(fig12[1, 1], aspect = AxisAspect(1))
ax4 = Axis(fig12[1, 1], aspect = AxisAspect(1))
ax3 = Axis(fig12[1, 1], aspect = AxisAspect(1))
for (i, ssp) in enumerate(ssp2100)
    vlines!(ax4, ssp, label = ssp_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssp_line_opts...)
    vlines!(ax4, ssp2100_hi[i], label = ssphi_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssphi_line_opts...)
    # vlines!(ax1, ssp, color = xpcolors[ssp_labels[i]]; ssp_line_opts...)
end

sref = 400
f_retreat = aqef_intermediate.f[end][1:sref:end] ./ polar_amplification .+ f2020
V_retreat = aqef_intermediate.V_sle[end][1:sref:end]
cmap = cgrad(:darktest, range(0, stop=1, length=100))
crange = (1, 11)
lines!(
    ax3,
    f_retreat,
    V_retreat,
    linewidth = lws[end],
    # label = xp_labels[end],
    color = f_retreat,
    colormap = cmap,
    colorrange = crange,
)
for i in 2:aqef_intermediate.n_xps-1
    f = aqef_intermediate.f[i][1:s:end] ./ polar_amplification .+ f2020
    V = aqef_intermediate.V_sle[i][1:s:end]
    p = (maximum(f) - crange[1]) / (crange[2] - crange[1])
    c = cmap[round(Int, p * 99) + 1]
    lines!(
        ax3,
        f,
        V,
        linewidth = lws[i],
        color = c,
        label = L"$f_\mathrm{max} = %$(round(f[1], digits = 1)) \, \mathrm{K}$",
    )
end
axislegend(ax3, position = :rt)
Legend(fig12[0, 1], ax4, nbanks = 4)
rowsize!(fig12.layout, 0, 50)
colsize!(fig12.layout, 1, 700)

format_axs!(ax3, atm_ax3, x1, x2)
hidedecorations!(atm_ax3)
hidedecorations!(ax4)
xlims!(ax4, 0, 11)
fig12
save(plotsdir("v2/hysteresis/fig12.png"), fig12)
save(plotsdir("v2/hysteresis/fig12.pdf"), fig12)