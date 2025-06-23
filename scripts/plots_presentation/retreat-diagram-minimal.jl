include("../intro.jl")

T = Float32
visc_type = "lowvisc"    # "equil" or "aqef"
heatmap_frames = "aqef"    # "equil" or "aqef"

xps = [
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-ocnforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-normvisc-fastatmforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-normforcing-withrestarts"),
]
xp_labels = [
    "OCN",
    "ATM",
    "REF",
]
aqef = AQEFResults(T, xps)
lws = [3, 3, 5]
cycling_colors = [:royalblue, :orange, :black]

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
xp_idx = aqef.n_xps
f2015 = 1.2

set_theme!(theme_latexfonts())
fig = Figure(size=(500, 400), fontsize = 18)
axs = [Axis(fig[i, j]) for i in 1:1, j in 1:1]
s = 400

shading = true

if shading
    large_misi_color = :red
    small_misi_color = :gray60
    shapeff_color = :gray80
    large_misi_alpha = 0.1
    small_misi_alpa = 0.2
    shapeff_alpha = 0.2

    shade_large_misi = [(1.2, 1.3), (4.7, 4.8), (6.1, 6.2), (7.0, 7.1), (7.8, 7.9)]
    for i in eachindex(shade_large_misi)
        if i == 1
            vlines!(axs[1, 1], (shade_large_misi[i][1]:0.01:shade_large_misi[i][2]) .+ f2015,
                alpha = large_misi_alpha, color = large_misi_color, label = "Large MISI")
        else
            vlines!(axs[1, 1], (shade_large_misi[i][1]:0.01:shade_large_misi[i][2]) .+ f2015,
                alpha = large_misi_alpha, color = large_misi_color)
        end
    end

    shade_small_misi = [(4.6, 4.7), (8.5, 8.6)]
    for i in eachindex(shade_small_misi)
        if i == 1
            vlines!(axs[1, 1], (shade_small_misi[i][1]:0.01:shade_small_misi[i][2]),
                alpha = small_misi_alpa, color = small_misi_color, label = "Small MISI")
        else
            vlines!(axs[1, 1], (shade_small_misi[i][1]:0.01:shade_small_misi[i][2]),
                alpha = small_misi_alpa, color = small_misi_color)
        end
    end

    shade_shapeff = [(9.9, 10), (10.2, 10.3)]
    for i in eachindex(shade_shapeff)
        if i == 1
            vlines!(axs[1, 1], (shade_shapeff[i][1]:0.01:shade_shapeff[i][2]),
                alpha = shapeff_alpha, color = shapeff_color, label = "Shape efficiency")
        else
            vlines!(axs[1, 1], (shade_shapeff[i][1]:0.01:shade_shapeff[i][2]),
                alpha = shapeff_alpha, color = shapeff_color)
        end
    end
end

if shading
    exp2plot = [aqef.n_xps]
else
    exp2plot = 1:aqef.n_xps
end

for k in exp2plot
    lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification .+ f2015, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k], color = color = cycling_colors[k])
end

axs[1, 1].xticks = 0:2:12
axs[1, 1].xminorticks = 0:0.2:12
axs[1, 1].yticks = 0:10:60
axs[1, 1].yminorticks = IntervalsBetween(10)
axs[1, 1].xminorgridvisible = true
axs[1, 1].yminorgridvisible = true
axs[1, 1].xaxisposition = :bottom
axs[1, 1].xlabel = L"GMT anomaly $f$ (K)"
axs[1, 1].ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(axs[1, 1], 0, 60)
xlims!(axs[1, 1], 0, 12)

axislegend(axs[1, 1], position = :lb, nbanks = 1)
save(plotsdir("16km/hysteresis/retreat-diagram-min-$visc_type-shading-$shading.png"), fig)