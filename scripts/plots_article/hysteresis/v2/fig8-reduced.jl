include("../../../intro.jl")

T = Float32
polar_amplification = 1.8
f2020 = 1.2
f_to = 0.25

xp = datadir("output/ais/v2/hyster/regrowth/aqef/refnomslow/0")
file1D = joinpath(xp, "yelmo1D.nc")
file2Dsm = joinpath(xp, "yelmo2Dsm.nc")
file2D = joinpath(xp, "yelmo2D.nc")

t1D = ncread(file1D, "time")
t2Dsm = ncread(file2Dsm, "time")
t2D = ncread(file2D, "time")
f = ncread(file1D, "hyst_f_now") ./ polar_amplification .+ f2020
x = ncread(file2Dsm, "xc")
y = ncread(file2Dsm, "yc")
# f_grnd = ncread(file2D, "f_grnd")
X, Y = ndgrid(x, y)

f_bif_rsb = 3.2
f_bif_asb = 2.7
n_grz = 20       # number of plotted grounding zones
di_grz = 2
dt_2D = mean(diff(t2D))
dt_grz = di_grz .* dt_2D
f_bif = f_bif_asb
i_bif = findlast(f .>= f_bif)
t_bif = t1D[i_bif]
t_end = t_bif .+ dt_grz .* n_grz
i_2D = argmin(abs.(t_bif .- t2D))
grz_steps = range(i_2D, step = di_grz, length = n_grz)

xl = (1400, 2500)
yl = (-1300, -200)

mask(X, Y, xl, yl) = (xl[1] .<= X .<= xl[2]) .& (yl[1] .<= Y .<= yl[2])
ase_mask = mask(X, Y, xl, yl)
mb = masked_massbalance(file2Dsm, t_bif - dt_2D, t_end, ase_mask)

# lines!(ax_bmb, mb.time ./ 1f3, mb.bmb, linewidth = lw2, label = "basal")

tot_grline_bmb = fill(NaN, n_grz)
tot_ocn_bmb = fill(NaN, n_grz)
tot_calving = fill(NaN, n_grz)
cnt_grline = fill(NaN, n_grz)
cnt_ocnmlt = fill(NaN, n_grz)
cnt_calving = fill(NaN, n_grz)

for l in eachindex(grz_steps)
    grline_mask = (0.5 .< ncslice(file2D, "mask_ocn", grz_steps[l]) .< 1.5) .& ase_mask
    ocnmlt_mask = (ncslice(file2D, "f_grnd", grz_steps[l]) .< 1) .&
        (ncslice(file2D, "f_ice", grz_steps[l]) .> 0.5) .&
        ase_mask
    calving_mask = (ncslice(file2D, "cmb", grz_steps[l]) .< 0) .& ase_mask
    cnt_grline[l] = count( grline_mask )
    cnt_ocnmlt[l] = count( ocnmlt_mask )
    cnt_calving[l] = count( calving_mask )
    tot_grline_bmb[l] = sum(ncslice(file2D, "bmb", grz_steps[l]) .* grline_mask)
    tot_ocn_bmb[l] = sum(ncslice(file2D, "bmb", grz_steps[l]) .* ocnmlt_mask)
    tot_calving[l] = sum(ncslice(file2D, "cmb", grz_steps[l]) .* calving_mask)
    @show cnt_grline[l] cnt_ocnmlt[l] cnt_calving[l]
end

mean_grline_bmb = tot_grline_bmb ./ cnt_grline
mean_ocn_bmb = tot_ocn_bmb ./ cnt_ocnmlt
mean_calving = tot_calving ./ cnt_calving

set_theme!(theme_latexfonts())
fs = 32
lw1, lw2 = 3, 6
fig8 = Figure(size = (1550, 1050), fontsize = fs)
ga = fig8[1, 1] = GridLayout()
ax_bmb = Axis(ga[1, 1])
ax_mbmb = Axis(ga[2, 1])
ax_cnt = Axis(ga[3, 1])
ylims!(ax_bmb, (-300, -30))
ylims!(ax_mbmb, (-1.5, -0.5))
ylims!(ax_cnt, (30, 500))
ax_bmb.yscale = Makie.pseudolog10
ax_cnt.yscale = Makie.pseudolog10
ax_hm = Axis(fig8[1, 2:3], aspect = DataAspect())
inset_ax = Axis(fig8[1, 2:3], width = Relative(0.3), height = Relative(0.3),
    halign = 0.98, valign = 0.98)
hidedecorations!(inset_ax)
heatmap!(inset_ax, x, y, ncslice(file2D, "z_bed", i_2D); cmaps["z_bed6"]...)
contour!(inset_ax, x, y,
    (xl[1] .< X .< xl[2]) .& (yl[1] .< Y .< yl[2]),
    levels = [0.5], color = :red, linewidth = 5)

catjet = cgrad(:jet, range(0, stop = 1, length = n_grz+1), categorical = true)
tsjet = cgrad(:jet)
hidedecorations!(ax_hm)
t = t1D[i_bif]
k = argmin(abs.(t .- t2D))
heatmap!(ax_hm, x, y, ncslice(file2D, "z_bed", k); cmaps["z_bed6"]...)
    grz_steps = range(k, step = di_grz, length = n_grz)
for ax in [ax_bmb, ax_mbmb, ax_cnt]
    hlines!(ax, 0, color = :black, linewidth = lw1, linestyle = :dash)
    vlines!(ax, t2D[grz_steps][1:2:end] ./ 1f3, color = [(c, 1) for c in catjet],
        linewidth = 3, linestyle = :dash)
end
for l in eachindex(grz_steps)[1:2:end]
    contour!(ax_hm, x, y,
        ncslice(file2D, "f_grnd", grz_steps[l]) .*
        (ncslice(file2D, "H_ice", grz_steps[l]) .> 100),
        levels = [0.5], color = catjet[l], linewidth = lw2)
end

# lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_ocn_bmb, linewidth = lw2)
# lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_calving, linewidth = lw2, color = :gray60)
lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_grline_bmb, linewidth = lw2)
# lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_ocn_bmb, linewidth = lw2, label = "basal, shelf")
# lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_calving, linewidth = lw2, label = "calving", color = :gray60)
lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_grline_bmb, linewidth = lw2, label = "grounding zone")
# lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_ocnmlt, linewidth = lw2)
# lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_calving, linewidth = lw2, color = :gray60)
lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_grline, linewidth = lw2)
# lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_calving, linewidth = lw2)
# lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_calving, linewidth = lw2, label = "calving")
# lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_calving, linewidth = lw2)

ax_bmb.ylabel = L"Total BMB ($\mathrm{m \, yr^{-1}}$)"
ax_mbmb.ylabel = L"Mean MB ($\mathrm{m \, yr^{-1}}$)"
ax_cnt.ylabel = "Cell count (1)"
ax_bmb.xticklabelsvisible = false
ax_mbmb.xticklabelsvisible = false
ax_cnt.xlabel = "Time (kyr)"

# heatmap!(ax_hm, x, y, ase_mask, colormap = cgrad([tr_color, tr_color]),
#     colorrange = (1, 1.1), lowclip = :transparent, highclip = :transparent)
xlims!(ax_hm, xl)
ylims!(ax_hm, yl)

ax_cnt.xlabel = "Time (kyr)"
# axislegend(ax_mbmb, nbanks = 1, fontsize = fs, position = :rc)
Colorbar(fig8[2, 2:3], label = "Bed elevation (km)", vertical = false,
    width = Relative(0.5), height = Relative(1), ticks = latexifyticks(-1:1, 1e3), halign = :right,
    flipaxis = false; cmaps["z_bed6"]...)
# Legend(fig8[2, 2:3], ax_mbmb, nbanks = 2, halign = 0.05, valign = 0.95)

ax_bmb.yticks = [-30, -45, -70, -100, -150, -230, -340][1:2:end]
ax_cnt.yticks = [30, 45, 70, 100, 150, 230, 340][1:2:end]
ax_bmb.yminorticks = IntervalsBetween(2)
ax_mbmb.yminorticks = IntervalsBetween(2)
ax_cnt.yminorticks = IntervalsBetween(2)
ax_bmb.yminorticksvisible = true
ax_mbmb.yminorticksvisible = true
ax_cnt.yminorticksvisible = true
dc = 450
colgap!(fig8.layout, 5)
rowgap!(fig8.layout, 5)
rowgap!(fig8.layout, 1, -60)
rowsize!(fig8.layout, 1, dc*2)
colsize!(fig8.layout, 1, dc)
colsize!(fig8.layout, 2, dc)
colsize!(fig8.layout, 3, dc)

xlims!(ax_bmb, (568, 587))
xlims!(ax_mbmb, (568, 587))
xlims!(ax_cnt, (568, 587))
text!(ax_bmb, 568.4, -30, text = "(a)", font = :bold)
text!(ax_mbmb, 568.4, -3, text = "(b)", font = :bold)
text!(ax_cnt, 568.4, 3, text = "(c)", font = :bold)
text!(ax_hm, xl[1]+20, yl[2]-50, text = "(d)", font = :bold)

fig8
save(plotsdir("v2/hysteresis/fig8-reduced.png"), fig8)
save(plotsdir("v2/hysteresis/fig8-reduced.pdf"), fig8)