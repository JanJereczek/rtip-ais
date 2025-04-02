include("../intro.jl")

T = Float32
visc_type = "normvisc"    # "equil" or "aqef"
xps = [
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-garbeforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-normforcing-withrestarts"),
    datadir("output/ais/hyster/16km/regrowth/aqef/pmpt-$visc_type-fastnormforcing-restarted"),
    datadir("output/ais/hyster/16km/regrowth/aqef/pmpt-$visc_type-fastnormforcing"),
]

xp_labels = [
    "HOW",
    "REF",
    nothing,
    nothing,
]
lws = [3, 5, 5, 5]
cycling_colors = [3, 1, 1, 1]
aqef = AQEFResults(T, xps)

polar_amplification = 1.8
f_to = 0.25
f_ref = aqef.f[end] ./ polar_amplification
stiching_forcing = 2
stiching_idx = findfirst(f_ref .<= stiching_forcing)

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 15
nrows, ncols = 2, 1
fig = Figure(size=(500, 1000), fontsize = 24)
axs = [Axis(fig[i, j], aspect = AxisAspect(1)) for i in 1:nrows, j in 1:ncols]
s = 50
for k in 1:aqef.n_xps
    if k < aqef.n_xps
        lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification, aqef.V_sle[k][1:s:end],
            linewidth = lws[k], label = xp_labels[k], color = Cycled(cycling_colors[k]))
    else
        lines!(axs[1, 1], aqef.f[k][1:s:stiching_idx] ./ polar_amplification,
            aqef.V_sle[k][1:s:stiching_idx], linewidth = lws[k], label = xp_labels[k],
            color = Cycled(cycling_colors[k]))
    end
end
# scatterlines!(axs[1, 1], eql.f ./ polar_amplification, eql.V_sle;
#     linewidth = lws[xp_idx], color = :black, label = "EQL", markersize = ms1)

g20dir = datadir("processed/garbe2020")
g20_retreat_ramp, _ = readdlm("$g20dir/retreat-ramp.csv", ',', header = true)
g20_retreat_equil, _ = readdlm("$g20dir/retreat-equil.csv", ',', header = true)
g20_regrowth_ramp, _ = readdlm("$g20dir/regrowth-ramp.csv", ',', header = true)
g20_regrowth_equil, _ = readdlm("$g20dir/regrowth-equil.csv", ',', header = true)
g20 = [g20_retreat_ramp, g20_retreat_equil, g20_regrowth_ramp, g20_regrowth_equil]
g20_labels = ["G20R", "G20E", nothing, nothing]
g20_colors = [5, 6, 5, 6]
g20_plotstyle = [:lines, :scatterlines, :lines, :scatterlines]
for i in eachindex(g20)
    if g20_plotstyle[i] == :lines
        lines!(axs[1, 1], g20[i][:, 1] ./ polar_amplification, g20[i][:, 2], linewidth = 3,
            label = g20_labels[i], color = Cycled(g20_colors[i]))
    elseif g20_plotstyle[i] == :scatterlines
        scatterlines!(axs[1, 1], g20[i][:, 1] ./ polar_amplification, g20[i][:, 2], linewidth = 3,
            label = g20_labels[i], color = Cycled(g20_colors[i]))
    end
end

discrete_support = collect(-5.0:15.0)
continuous_support = collect(range(-5, stop = 15, step = 0.1))

function keep_unique(x, y)
    i = unique(i -> x[i], eachindex(x))
    x = x[i]
    y = y[i]
    return x, y
end

function diff_on_support(x1, y1, x2, y2, support)

    x1, y1 = keep_unique(x1, y1)
    x2, y2 = keep_unique(x2, y2)

    idx1 = sortperm(x1)
    idx2 = sortperm(x2)
    y1_interp = linear_interpolation(x1[idx1], y1[idx1], extrapolation_bc = NaN)
    y2_interp = linear_interpolation(x2[idx2], y2[idx2], extrapolation_bc = NaN)
    Y1 = y1_interp.(support)
    Y2 = y2_interp.(support)
    dy = Y1 .- Y2
    return Y1, Y2, dy
end

y1ref, y2ref, dref = diff_on_support(aqef.f[2] ./ polar_amplification, aqef.V_sle[2],
    aqef.f[3] ./ polar_amplification, aqef.V_sle[3], continuous_support)
# dhow = diff_on_support(aqef.f[1], aqef.V_sle[1], aqef.f[1], aqef.V_sle[1], continuous_support)
y1g20r, y2g20r, dg20r = diff_on_support(g20_retreat_ramp[:, 1] ./ polar_amplification,
    g20_retreat_ramp[:, 2], g20_regrowth_ramp[:, 1] ./ polar_amplification,
    g20_regrowth_ramp[:, 2], continuous_support)
y1g20e, y2g20e, dg20e = diff_on_support(g20_retreat_equil[:, 1] ./ polar_amplification,
    g20_retreat_equil[:, 2], g20_regrowth_equil[:, 1] ./ polar_amplification,
    g20_regrowth_equil[:, 2], discrete_support)

lines!(axs[2, 1], continuous_support, dref, linewidth = 3, label = "REF", color = Cycled(1))
# lines!(axs[2, 1], continuous_support, dhow, linewidth = 3, label = "HOW", color = Cycled(3))
lines!(axs[2, 1], continuous_support, dg20r, linewidth = 3, label = "G20R", color = Cycled(5))
scatterlines!(axs[2, 1], discrete_support, dg20e, linewidth = 3, label = "G20E", color = Cycled(6))

# lines!(axs[1, 1], continuous_support, y1ref, linestyle = :dash, linewidth = 3, color = Cycled(2))
# lines!(axs[1, 1], continuous_support, y2ref, linestyle = :dash, linewidth = 3, color = Cycled(2))
# # lines!(axs[1, 1], continuous_support, y1how, linestyle = :dash, linewidth = 3, color = Cycled(3))
# # lines!(axs[1, 1], continuous_support, y2how, linestyle = :dash, linewidth = 3, color = Cycled(3))
# lines!(axs[1, 1], continuous_support, y1g20r, linestyle = :dash, linewidth = 3, color = Cycled(3))
# lines!(axs[1, 1], continuous_support, y2g20r, linestyle = :dash, linewidth = 3, color = Cycled(3))
# lines!(axs[1, 1], discrete_support, y1g20e, linestyle = :dash, linewidth = 3, color = Cycled(4))
# lines!(axs[1, 1], discrete_support, y2g20e, linestyle = :dash, linewidth = 3, color = Cycled(4))

axs[1, 1].xticks = -10:2:20
axs[2, 1].xticks = -10:2:20
axs[1, 1].xminorticks = -10:0.2:20
axs[2, 1].xminorticks = -10:0.2:20
axs[1, 1].xminorgridvisible = true
axs[2, 1].xminorgridvisible = true
axs[1, 1].xaxisposition = :top
axs[1, 1].xlabel = L"GMT anomaly $f$ (K)"
axs[2, 1].xlabel = L"GMT anomaly $f$ (K)"

axs[1, 1].yticks = 0:10:60
axs[1, 1].yminorticks = IntervalsBetween(10)
axs[1, 1].yminorgridvisible = true
axs[1, 1].ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
axs[2, 1].ylabel = L"$\Delta V_\mathrm{af}$ (m SLE)"
axislegend(axs[1, 1], position = :lb)
axislegend(axs[2, 1], position = :lt)
ylims!(axs[1, 1], 0, 60)
xlims!(axs[1, 1], -5, 12)
xlims!(axs[2, 1], -5, 12)

save(plotsdir("16km/hysteresis/retreat-regrowth.png"), fig)