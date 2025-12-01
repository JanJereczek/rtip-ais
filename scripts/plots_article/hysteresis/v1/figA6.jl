include("../../intro.jl")

T = Float32
polar_amplification = 1.8
f_to = 0.25

xp = datadir("output/ais/hyster/16km/retreat/aqef/pmpt-lowvisc-normforcing-withrestarts/0")
file1D = joinpath(xp, "yelmo1D.nc")
file2Dsm = joinpath(xp, "yelmo2Dsm.nc")
file2D = joinpath(xp, "yelmo2D.nc")

t1D = ncread(file1D, "time")
t2Dsm = ncread(file2Dsm, "time")
t2D = ncread(file2D, "time")
f = ncread(file1D, "hyst_f_now") ./ polar_amplification
x = ncread(file2Dsm, "xc")
y = ncread(file2Dsm, "yc")
# f_grnd = ncread(file2D, "f_grnd")
X, Y = ndgrid(x, y)

n_grz = 12       # number of plotted grounding zones
di_grz = 1
dt_2D = mean(diff(t2D))
dt_grz = di_grz .* dt_2D
f_bif = 1.2
i_bif = findlast(f .<= f_bif) .+ 2
t_bif = t1D[i_bif]
t_end = t_bif .+ dt_grz .* n_grz

xl = (-2000, 0)
yl = (-1400, 600)

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
fs = 38
lw1, lw2 = 3, 6
fig = Figure(size = (1550, 1050), fontsize = fs)
ga = fig[1, 1] = GridLayout()
ax_bmb = Axis(ga[1, 1])
ax_mbmb = Axis(ga[2, 1])
ax_cnt = Axis(ga[3, 1])
ax_hm = Axis(fig[1, 2:3], aspect = DataAspect())

catjet = cgrad(:jet, range(0, stop = 1, length = n_grz+1), categorical = true)
tsjet = cgrad(:jet)
hidedecorations!(ax_hm)
t = t1D[i_bif]
k = argmin(abs.(t .- t2D))
heatmap!(ax_hm, x, y, ncslice(file2D, "z_bed", k); cmaps["z_bed6"]...)
    grz_steps = range(k, step = di_grz, length = n_grz)
for ax in [ax_bmb, ax_mbmb, ax_cnt]
    hlines!(ax, 0, color = :black, linewidth = lw1, linestyle = :dash)
    vlines!(ax, t2D[grz_steps] ./ 1f3, color = [(c, 0.4) for c in catjet],
        linewidth = 6)
end
for l in eachindex(grz_steps)
    contour!(ax_hm, x, y,
        ncslice(file2D, "f_grnd", grz_steps[l]) .*
        (ncslice(file2D, "H_ice", grz_steps[l]) .> 100),
        levels = [0.5], color = catjet[l], linewidth = lw2)
end

lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_ocn_bmb, linewidth = lw2)
lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_grline_bmb, linewidth = lw2)
lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_ocn_bmb, linewidth = lw2, label = "shelf")
lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_grline_bmb, linewidth = lw2, label = "gr. zone")
lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_ocnmlt, linewidth = lw2)
lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_grline, linewidth = lw2)
# lines!(ax_bmb, t2D[grz_steps] ./ 1f3, tot_calving, linewidth = lw2)
# lines!(ax_mbmb, t2D[grz_steps] ./ 1f3, mean_calving, linewidth = lw2, label = "calving")
# lines!(ax_cnt, t2D[grz_steps] ./ 1f3, cnt_calving, linewidth = lw2)

ax_bmb.ylabel = L"Tot. BMB ($\mathrm{m \, yr^{-1}}$)"
ax_mbmb.ylabel = L"Mean BMB ($\mathrm{m \, yr^{-1}}$)"
ax_cnt.ylabel = "Cell count (1)"
ax_bmb.xticklabelsvisible = false
ax_mbmb.xticklabelsvisible = false
ax_cnt.xlabel = "Time (kyr)"

# heatmap!(ax_hm, x, y, ase_mask, colormap = cgrad([tr_color, tr_color]),
#     colorrange = (1, 1.1), lowclip = :transparent, highclip = :transparent)
xlims!(ax_hm, xl)
ylims!(ax_hm, yl)

ax_cnt.xlabel = "Time (kyr)"
axislegend(ax_mbmb, nbanks = 1, fontsize = fs, position = :rb)
Colorbar(fig[2, 2:3], label = "Bed elevation (km)", vertical = false,
    width = Relative(0.5), height = Relative(1), ticks = latexifyticks(-1:1, 1e3),
    flipaxis = false; cmaps["z_bed6"]...)

dc = 450
colgap!(fig.layout, 5)
rowgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -80)
rowsize!(fig.layout, 1, dc*2)
colsize!(fig.layout, 1, dc)
colsize!(fig.layout, 2, dc)
colsize!(fig.layout, 3, dc)
save(plotsdir("16km/hysteresis/figA6.png"), fig)
save(plotsdir("16km/hysteresis/figA6.pdf"), fig)