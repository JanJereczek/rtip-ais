include("../../../intro.jl")

sr = StepRtip{Float32}(visc_cases)
prefix = "output/ais/v2/ramps/"
dirs = [
    datadir("$prefix/steps-sigmarange"),
]
for dir in dirs
    aggregate_xp!(sr, dir)
end
sort!(sr)

basins = ["WAIS", "WSB", "RSB", "ASB ext.", "ASB int."]
f_gmt_bif = [1.97, 5.89, 7.0, 7.97, 8.67]
V_bif = [54, 42, 31, 20, 13.5]
i_rtip = [[findlast(V .> V_bif[i]) + 1 for i in eachindex(f_gmt_bif)] for V in sr.V]
f_rtip = [sr.f[i][i_rtip[i]] for i in eachindex(i_rtip)]

# f_rtip = Float32[1.3003739, 8.431313, 10.136125, 11.998654, 13.31169]
# f_rtip = Float32[1.3774256, 8.441516, 10.356382, 12.117172, 13.441652]
# f_rtip = Float32[1.4284503, 8.41603, 10.551467, 12.18136, 13.470692]

for i in eachindex(f_rtip)
    f_rtip[i][2] = max.(f_rtip[i][2], 5.75)
end

eqfig = Figure(size=(1600, 900), fontsize = 20)
axeq = Axis(eqfig[1, 1])
for i in eachindex(sr.visc_cases)
    lines!(
        axeq,
        sr.f[i],
        sr.V[i],
        label = visc_labels[i];
        linewidth = 4,
        color = viscmap[i],
    )
    vlines!(axeq, f_rtip[i], color = viscmap[i], linestyle = :dash)
end
axislegend(axeq, position = :rt)
vlines!(axeq, f_gmt_bif, color = :black, linestyle = :dash, linewidth = 4,
    label = "Bifurcation points")
save(plotsdir("v2/rtip/equil.png"), eqfig)
save(plotsdir("v2/rtip/equil.pdf"), eqfig)

lw = 4
ms = 30
slines_opts = (linewidth = lw,)
set_theme!(theme_latexfonts())
fig = Figure(size=(1000, 600), fontsize = 20)

ax = Axis(fig[1, 1], aspect = AxisAspect(1))
tbl_bif = (
    cat = repeat(1:5, outer = 5),
    height = vcat(1.1 * f_rtip...),
    grp = repeat(1:5, inner = 5),
)
tbl_rtip = (
    cat = repeat(1:5, outer = 5),
    height = vcat(f_rtip...),
    grp = repeat(1:5, inner = 5),
)
barplot!(
    ax,
    tbl_bif.cat,
    tbl_bif.height,
    dodge = tbl_bif.grp,
    color = viscmap[tbl_bif.grp],
    alpha = 0.5,
)
barplot!(
    ax,
    tbl_rtip.cat,
    tbl_rtip.height,
    dodge = tbl_rtip.grp,
    color = viscmap[tbl_rtip.grp],
    alpha = 1,
)

ax.xticks = (eachindex(basins), basins)
ax.xticklabelrotation = pi/8
ax.ylabel = L"$f^\mathrm{crit}$ (K)"
ax.yticks = 0:2:10
ylims!(ax, (1.5, 9))

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
for (i, ssp) in enumerate(ssp2100)
    hlines!(ax, ssp, label = ssp_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssp_line_opts...)
end
for (i, ssphi) in enumerate(ssp2100_hi)
    hlines!(ax, ssphi, label = ssphi_short_labels[i], color = xpcolors[ssp_labels[i]];
        ssphi_line_opts...)
end
fig

ax2 = Axis(fig[1, 2], aspect = AxisAspect(1))
hlines!(ax2, 0, color = :gray70, linestyle = :dash, linewidth = lw)
tbl_gap = (
    cat = repeat(1:5, outer = 5),
    height = vcat([f_rtip[i] .- f_gmt_bif for i in eachindex(f_rtip)]...),
    grp = repeat(1:5, inner = 5),
)

barplot!(ax2, tbl_gap.cat, tbl_gap.height,
        dodge = tbl_gap.grp,
        color = viscmap[tbl_gap.grp])

ax2.xticks = (eachindex(basins), basins)
ax2.xticklabelrotation = pi/8
ax2.ylabel = L"$\gamma$ (K)"
ax2.yticks = -1.6:0.2:0.1
ax2.yaxisposition = :right
yl = (-0.5, 0.1)
ylims!(ax2, yl)

Legend(fig[0, 1], ax, nbanks = 4)
labels = vcat(visc_labels...)
elements = [PolyElement(polycolor = viscmap[i]) for i in 1:5]
Legend(fig[0, 2], elements, labels, nbanks = 5)

rowsize!(fig.layout, 0, 60)
rowgap!(fig.layout, 1, -10)
colsize!(fig.layout, 1, 420)
colsize!(fig.layout, 2, 420)
colgap!(fig.layout, 5)
fig



# save(plotsdir("v2/rtip/gap.png"), fig)
# save(plotsdir("v2/rtip/gap.pdf"), fig)