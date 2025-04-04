include("../intro.jl")

T = Float32
visc_type = "normvisc"    # "equil" or "aqef"
xps = [
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-garbeforcing"),
    datadir("output/ais/hyster/16km/regrowth/aqef/pmpt-$visc_type-garbeforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-fastatmforcing"),
    datadir("output/ais/hyster/16km/regrowth/aqef/pmpt-$visc_type-atmforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-normforcing-withrestarts"),
    datadir("output/ais/hyster/16km/regrowth/aqef/pmpt-$visc_type-fastnormforcing-restarted"),
    datadir("output/ais/hyster/16km/regrowth/aqef/pmpt-$visc_type-fastnormforcing"),
]

xp_labels = [
    "HOW",
    nothing,
    "ATM",
    nothing,
    "REF",
    nothing,
    nothing,
]
aqef = AQEFResults(T, xps)

eqldir1 = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
eql1 = EquilResults(T, eqldir1)
eqldir2 = datadir("output/ais/hyster/16km/regrowth/equil")
eql2 = EquilResults(T, eqldir2)

lw1, lw2, lw3 = 3, 5, 7
ms1, ms2 = 15, 25
lws = [lw2, lw2, lw2, lw2, lw3, lw3, lw3]
cycling_colors = [:midnightblue, :midnightblue, 2, 2, :steelblue1, :steelblue1, :steelblue1]
h94c = :gray80
g20rc = :gray60
g20ec = :gray30

polar_amplification = 1.8
f_to = 0.25
f_ref = aqef.f[end] ./ polar_amplification
stiching_forcing = 2
stiching_idx = findfirst(f_ref .<= stiching_forcing)

set_theme!(theme_latexfonts())
fig = Figure(size=(700, 700), fontsize = 24)
ssp_ax = Axis(fig[1, 1], aspect = AxisAspect(1))
hidedecorations!(ssp_ax)
ax = Axis(fig[1, 1], aspect = AxisAspect(1))
atm_ax = Axis(fig[1, 1], aspect = AxisAspect(1))
hideydecorations!(atm_ax)

hist = readdlm(datadir("processed/SSP/History.csv"), ',')
f2014 = hist[end, 2]
ssp1 = readdlm(datadir("processed/SSP/SSP1.csv"), ',') # .- f2014
ssp2 = readdlm(datadir("processed/SSP/SSP2.csv"), ',') # .- f2014
ssp3 = readdlm(datadir("processed/SSP/SSP3.csv"), ',') # .- f2014
ssp5 = readdlm(datadir("processed/SSP/SSP5.csv"), ',') # .- f2014
ssp1_2100 = ssp1[end, 2]
ssp2_2100 = ssp2[end, 2]
ssp3_2100 = ssp3[end, 2]
ssp5_2100 = ssp5[end, 2]

line_opts = (linewidth = lw1, linestyle = :dash)
vlines!(ssp_ax, [ssp1_2100], color = :darkblue, label = "SSP1-2100"; line_opts...)
vlines!(ssp_ax, [ssp2_2100], color = :lightblue, label = "SSP2-2100"; line_opts...)
vlines!(ssp_ax, [ssp3_2100], color = :orange, label = "SSP3-2100"; line_opts...)
vlines!(ssp_ax, [ssp5_2100], color = :darkred, label = "SSP5-2100"; line_opts...)

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
g20_colors = [g20rc, g20ec, g20rc, g20ec]
g20_plotstyle = [:lines, :scatterlines, :lines, :scatterlines]
g20_markers = [nothing, :dtriangle, nothing, :utriangle]
for i in eachindex(g20)
    if g20_plotstyle[i] == :lines
        lines!(ax, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2], linewidth = 3,
            label = g20_labels[i], color = g20_colors[i])
    elseif g20_plotstyle[i] == :scatterlines
        scatterlines!(ax, g20[i][:, 1] ./ polar_amplification, g20[i][:, 2],
            linewidth = 3, label = g20_labels[i], color = g20_colors[i],
            marker = g20_markers[i], markersize = ms1)
    end
end

m3_msle = 1e9 * 1e6 / (3.625 * 1e14)
# h94[1, 2] * m3_msle * gr_correction = 58
gr_correction = 58 / (h94_retreat[1, 2] * m3_msle)

scatterlines!(ax, h94_retreat[:, 1] ./ polar_amplification,
    h94_retreat[:, 2] .* m3_msle .* gr_correction, marker = :dtriangle,
    linewidth = lw1, label = "H94", color = h94c, markersize = ms1)
scatterlines!(ax, h94_regrowth[:, 1] ./ polar_amplification,
    h94_regrowth[:, 2] .* m3_msle .* gr_correction, marker = :utriangle,
    linewidth = lw1, color = h94c, markersize = ms1)


s = 400
for k in 1:aqef.n_xps
    if k < aqef.n_xps
        lines!(ax, aqef.f[k][1:s:end] ./ polar_amplification .+ f2014,
            aqef.V_sle[k][1:s:end], linewidth = lws[k], label = xp_labels[k],
            color = lcolor(cycling_colors[k]))
    else
        lines!(ax, aqef.f[k][1:s:stiching_idx] ./ polar_amplification .+ f2014,
            aqef.V_sle[k][1:s:stiching_idx], linewidth = lws[k], label = xp_labels[k],
            color = lcolor(cycling_colors[k]))
    end
end
scatter!(ax, eql1.f ./ polar_amplification .+ f2014, eql1.V_sle;
    color = :black, label = "EQL", markersize = ms1)
scatter!(ax, eql2.f ./ polar_amplification .+ f2014, eql2.V_sle;
    color = :black, markersize = ms1)

ax.xticks = -10:2:20
ax.xminorticks = -10:0.2:20
ax.xminorgridvisible = true
ax.xlabel = L"GMT anomaly $f$ (K)"
ax.yticks = 0:10:60
ax.yminorticks = IntervalsBetween(10)
ax.yminorgridvisible = true
ax.ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"

atm_ax.xticks = -10:2:30
atm_ax.xgridvisible = false
atm_ax.xminorticks = -10:1:30
atm_ax.xaxisposition = :top
atm_ax.xlabel = L"Regional atmospheric temperature anomaly $f_a$ (K)"
axislegend(ax, position = :lb)
axislegend(ssp_ax, position = :rt)
ylims!(ax, 0, 60)
x1, x2 = 0, 12
xlims!(ax, x1, x2)
xlims!(ssp_ax, x1, x2)
xlims!(atm_ax, x1 * polar_amplification, x2 * polar_amplification)
save(plotsdir("16km/hysteresis/retreat-regrowth.png"), fig)


#=
axs[2, 1].ylabel = L"$\Delta V_\mathrm{af}$ (m SLE)"
axs[2, 1].xticks = -10:2:20
axs[2, 1].xminorticks = -10:0.2:20
axs[2, 1].xminorgridvisible = true
axs[2, 1].xlabel = L"GMT anomaly $f$ (K)"
axislegend(axs[2, 1], position = :lt)
xlims!(axs[2, 1], -5, 12)


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

# lines!(ax, continuous_support, y1ref, linestyle = :dash, linewidth = 3, color = Cycled(2))
# lines!(ax, continuous_support, y2ref, linestyle = :dash, linewidth = 3, color = Cycled(2))
# # lines!(ax, continuous_support, y1how, linestyle = :dash, linewidth = 3, color = Cycled(3))
# # lines!(ax, continuous_support, y2how, linestyle = :dash, linewidth = 3, color = Cycled(3))
# lines!(ax, continuous_support, y1g20r, linestyle = :dash, linewidth = 3, color = Cycled(3))
# lines!(ax, continuous_support, y2g20r, linestyle = :dash, linewidth = 3, color = Cycled(3))
# lines!(ax, discrete_support, y1g20e, linestyle = :dash, linewidth = 3, color = Cycled(4))
# lines!(ax, discrete_support, y2g20e, linestyle = :dash, linewidth = 3, color = Cycled(4))
=#