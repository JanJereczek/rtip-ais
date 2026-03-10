include("../../../intro.jl")

T = Float32
regrowth_dir = datadir("output/ais/v2/hyster/regrowth")
xps = [
    "$regrowth_dir/aqef/refnomslow-restarted",
    "$regrowth_dir/aqef/refnomslow",
]
aqef = AQEFResults(T, xps)
eql = EquilResults(T, "$regrowth_dir/equil/refnomslow")

xp_labels = [
    nothing,
    "REF",
]
lw1, lw2 = 3, 6
lws = [lw2, lw2]
cycling_colors = [
    xpcolors["REF"],
    xpcolors["REF"],
]

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
f2020 = 1.2
heatmap_frames = "aqef"    # "equil" or "aqef"
xp_idx = aqef.n_xps
f_ref = aqef.f[end] ./ polar_amplification
cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)

g20dir = datadir("processed/garbe2020")
g20_regrowth_ramp, _ = readdlm("$g20dir/regrowth-ramp.csv", ',', header = true)
g20_regrowth_equil, _ = readdlm("$g20dir/regrowth-equil.csv", ',', header = true)
g20_regrowth_equil = vcat([-3 55.6; -2 55.6; -1 55.3; 0 55], g20_regrowth_equil)

g20 = [g20_regrowth_ramp, g20_regrowth_equil]
g20_labels = ["G20R", "G20E"]
g20_colors = [xpcolors["G20R"], xpcolors["G20E"]]
g20_plotstyle = [:lines, :scatterlines]
g20_markers = [nothing, :utriangle]

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 18
fig4 = Figure(size=(800, 800), fontsize = 24)
ax = Axis(fig4[1, 1])
ax.xticks = -2:1:12
ax.xminorticks = -2:0.2:12
ax.yticks = 0:10:60
ax.yminorticks = IntervalsBetween(10)
ax.xminorgridvisible = true
ax.yminorgridvisible = true
ax.xlabel = L"GMT anomaly $f$ (K)"
ax.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(ax, 0, 60)
xlims!(ax, 0, 10)

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
    (1.174, "Wilkes", 35, 7.5),
    (0.9, "Wilkes & Siple Coast", 5, 22),
    (0.45, "Wilkes", 40, 7.5),
]

f_bif = [bifs[i][1] for i in eachindex(bifs)]
vlines!(ax, f_bif, color = :gray60, alpha = 0.5, linewidth = 5)

for i in eachindex(g20)
    if g20_plotstyle[i] == :lines
        lines!(ax, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2], linewidth = 3,
            label = g20_labels[i], color = g20_colors[i])
    elseif g20_plotstyle[i] == :scatterlines
        scatterlines!(ax, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2],
            linewidth = 3, label = g20_labels[i], color = g20_colors[i],
            marker = g20_markers[i], markersize = ms2)
    end
end
for i in eachindex(bifs)
    rectangle!(ax, (bifs[i][1] - 0.13, bifs[i][3] - 0.4), 0.3, 0.7*bifs[i][4])
    text!(ax, bifs[i][1] - 0.12, bifs[i][3], text = bifs[i][2], rotation = π/2, fontsize = 18)
end
vlines!(ax, 1f6, alpha = 0.9, color = :gray60, linewidth = 6, label = "Bifurcation")

s = 200
alpha = 1
for k in 1:aqef.n_xps
    lines!(ax, aqef.f[k][1:s:end] ./ polar_amplification .+ f2020,
        aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
        color = lcolor(cycling_colors[k]), alpha = alpha)
end

idx = sortperm(eql.f)
f = eql.f[idx]
V_sle = eql.V_sle[idx]
V_sle = filter_outliers(V_sle, mode = :low)
scatter!(ax, f ./ polar_amplification .+ f2020, V_sle;
    color = :cornflowerblue, label = "EQL", markersize = ms1)
axislegend(ax, position = :lt)
# Legend(fig4[0, 1], ax, nbanks = 5, framevisible = false)

# Compute mean diff between eql and ref
idx = sortperm(aqef.f[end])
aqef_itp = linear_interpolation(aqef.f[end][idx] ./ polar_amplification .+ f2020, aqef.V_sle[end][idx])
# aqef_eql = aqef_itp.(f ./ polar_amplification .+ f2020)
# mean_diff = mean(abs.(aqef_eql[f .< 16] .- V_sle[f .< 16])) # exclude ranges where volue = 0 to avoid underestimating error
# println("Mean difference REF vs EQL: $(round(mean_diff, digits=2)) m SLE")
# lines(f[f .< 16], aqef_eql[f .< 16] .- V_sle[f .< 16])
ax.xreversed = true

# rowsize!(fig4.layout, 0, 20)
colsize!(fig4.layout, 1, 700)
# axislegend(ax, position = :ct)
fig4
save(plotsdir("v2/hysteresis/fig4-reduced.png"), fig4)
save(plotsdir("v2/hysteresis/fig4-reduced.pdf"), fig4)