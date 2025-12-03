include("../../../intro.jl")

polar_amplification = 1.8
f_to = 0.25
f2015 = 1.2

dir = datadir("output/ais/v2/hyster")
file_retreat = "$dir/retreat/aqef/minvisc/refnomslow/0/yelmo2D.nc"
file_regrowth = "$dir/regrowth/aqef/refnomslow/0/yelmo2D.nc"
file_regrowth_1D = replace(file_regrowth, "yelmo2D.nc" => "yelmo1D.nc")
x, y = ncread(file_retreat, "xc"), ncread(file_retreat, "yc")
nx, ny = length(x), length(y)
f = ncread(file_regrowth_1D, "hyst_f_now") ./ polar_amplification .+ f2015
t_1D = ncread(file_regrowth_1D, "time")
t_2D = ncread(file_regrowth, "time")
f_intermediate = 4
k_intermediate_1D = argmin((f .- f_intermediate).^2)
k_intermediate = argmin((t_1D[k_intermediate_1D] .- t_2D).^2)

z_bed_t0 = ncread(file_regrowth, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])
z_bed_t1 = ncread(file_regrowth, "z_bed", start = [1, 1, k_intermediate], count = [-1, -1, 1])
z_bed_t2 = ncread(file_retreat, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])

xy_max = 2800
i1_ais = findfirst(x .> -xy_max)+10
i2_ais = findlast(x .< xy_max)
j1_ais = findfirst(y .> -xy_max)+20
j2_ais = findlast(y .< xy_max)-20

bsl_0 = 0
bsl_1 = 58
lw1 = 2
lw2 = 5

nrows = 1
ncols = 3
set_theme!(theme_latexfonts())
fig8 = Figure(size = (1500, 600), fontsize = 22)
axs = [Axis(fig8[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
axs[1, 1].title = "(a) Without ice (UPL)"
axs[1, 2].title = "(b) Intermediate regrowth"
axs[1, 3].title = "(c) Present-day ice thickness (DPR)"
z_bed = [z_bed_t0, z_bed_t1, z_bed_t2]
f_frames = [12, f_intermediate, f2015]
for i in eachindex(z_bed)
    zb = z_bed[i]
    heatmap!(axs[1, i], view(zb, i1_ais:i2_ais, j1_ais:j2_ais); cmaps["z_bed5"]...)
    hidedecorations!(axs[1, i])
    text!(axs[1, i], 5, 5, text="f = $(f_frames[i]) K", font = :bold, color = :white)
end

Colorbar(fig8[nrows+1, 1:3], label = "Bed elevation (km)", width = Relative(0.3),
    flipaxis = false, vertical = false, ticks = latexifyticks(-4:2, 1e3);
    cmaps["z_bed5"]...)

rowgap!(fig8.layout, -10)
colgap!(fig8.layout, 10)
fig8
save(plotsdir("v2/hysteresis/fig8.png"), fig8)
save(plotsdir("v2/hysteresis/fig8.pdf"), fig8)