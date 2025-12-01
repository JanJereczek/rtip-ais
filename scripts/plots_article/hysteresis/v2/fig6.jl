include("../../intro.jl")

file_t0 = datadir("output/ais/hyster/16km/retreat/aqef/"*
    "pmpt-normvisc-normforcing-withrestarts/0/yelmo2D.nc")
file_t2 = datadir("output/ais/hyster/16km/regrowth/aqef/"*
    "pmpt-normvisc-fastnormforcing/0/yelmo2D.nc")
x, y = ncread(file_t0, "xc"), ncread(file_t0, "yc")
nx, ny = length(x), length(y)
nt1 = length(ncread(file_t0, "time"))

z_bed_t0 = ncread(file_t2, "z_bed", start = [1, 1, 169], count = [-1, -1, 1])
z_bed_t1 = ncread(file_t0, "z_bed", start = [1, 1, nt1], count = [-1, -1, 1])
z_bed_t2 = ncread(file_t0, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])

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
fig = Figure(size = (1500, 600), fontsize = 22)
axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(ax) for ax in axs]
axs[1, 1].title = "(a) Intermediate regrowth"
axs[1, 2].title = "(b) Without ice (UPL)"
axs[1, 3].title = "(c) Present day (DPR)"
z_bed = [z_bed_t0, z_bed_t1, z_bed_t2]
for i in eachindex(z_bed)
    zb = z_bed[i]
    heatmap!(axs[1, i], view(zb, i1_ais:i2_ais, j1_ais:j2_ais); cmaps["z_bed5"]...)
end

Colorbar(fig[nrows+1, 1:3], label = "Bed elevation (km)", width = Relative(0.3),
    flipaxis = false, vertical = false, ticks = latexifyticks(-4:2, 1e3);
    cmaps["z_bed5"]...)

rowgap!(fig.layout, -10)
colgap!(fig.layout, 10)
save(plotsdir("16km/hysteresis/fig6.png"), fig)
save(plotsdir("16km/hysteresis/fig6.pdf"), fig)