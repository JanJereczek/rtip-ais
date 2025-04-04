include("../intro.jl")

file_t0 = datadir("output/ais/hyster/16km/retreat/aqef/"*
    "pmpt-normvisc-normforcing-withrestarts/0/yelmo2D.nc")
x, y = ncread(file_t0, "xc"), ncread(file_t0, "yc")
nx, ny = length(x), length(y)
nt1 = length(ncread(file_t0, "time"))
z_bed_t0 = ncread(file_t0, "z_bed", start = [1, 1, 1], count = [-1, -1, 1])
z_bed_t1 = ncread(file_t0, "z_bed", start = [1, 1, nt1], count = [-1, -1, 1])

xy_max = 2800
i1_ais = findfirst(x .> -xy_max)
i2_ais = findlast(x .< xy_max)
j1_ais = findfirst(y .> -xy_max)
j2_ais = findlast(y .< xy_max)

i1_wais = findfirst(x .> -2000)
i2_wais = findlast(x .< -400)
j1_wais = findfirst(y .> -1200)
j2_wais = findlast(y .< 400)

i1_asb = findfirst(x .> 200)
i2_asb = findlast(x .< 2700)
j1_asb = findfirst(y .> -2300)
j2_asb = findlast(y .< 200)

mask = zeros(eltype(z_bed_t0), size(z_bed_t0))
view(mask, i1_wais:i2_wais, j1_wais:j2_wais) .= 1
view(mask, i1_asb:i2_asb, j1_asb:j2_asb) .= 2

bsl_0 = 0
bsl_1 = 58
lw1 = 2
lw2 = 5

framecolors = [:red, :orange, :gray20, :gray70]
nrows = 2
ncols = 3
set_theme!(theme_latexfonts())
fig = Figure(size = (1500, 1120), fontsize = 32)
axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(ax) for ax in axs]
axs[1, 1].title = "Antarctica"
axs[1, 2].title = "West-Antarctica"
axs[1, 3].title = "Aurora and Wilkes"
axs[1, 1].ylabel = "Present day (DPR)"
axs[2, 1].ylabel = "Without ice (UPL)"
axs[1, 1].ylabelvisible = true
axs[2, 1].ylabelvisible = true
heatmap!(axs[1, 1], view(z_bed_t0, i1_ais:i2_ais, j1_ais:j2_ais); cmaps["z_bed5"]...)
heatmap!(axs[2, 1], view(z_bed_t1, i1_ais:i2_ais, j1_ais:j2_ais); cmaps["z_bed5"]...)
contour!(axs[1, 1], view(mask, i1_ais:i2_ais, j1_ais:j2_ais), levels = [0.5], linewidth = lw2,
    colormap = cgrad([framecolors[1], framecolors[1]]))
contour!(axs[1, 1], view(mask, i1_ais:i2_ais, j1_ais:j2_ais), levels = [1.5], linewidth = lw2,
    colormap = cgrad([framecolors[2], framecolors[2]]))
contour!(axs[2, 1], view(mask, i1_ais:i2_ais, j1_ais:j2_ais), levels = [0.5], linewidth = lw2,
    colormap = cgrad([framecolors[3], framecolors[3]]))
contour!(axs[2, 1], view(mask, i1_ais:i2_ais, j1_ais:j2_ais), levels = [1.5], linewidth = lw2,
    colormap = cgrad([framecolors[4], framecolors[4]]))
heatmap!(axs[1, 2], view(z_bed_t0, i1_wais:i2_wais, j1_wais:j2_wais); cmaps["z_bed5"]...)
heatmap!(axs[2, 2], view(z_bed_t1, i1_wais:i2_wais, j1_wais:j2_wais); cmaps["z_bed5"]...)
heatmap!(axs[1, 3], view(z_bed_t0, i1_asb:i2_asb, j1_asb:j2_asb); cmaps["z_bed5"]...)
heatmap!(axs[2, 3], view(z_bed_t1, i1_asb:i2_asb, j1_asb:j2_asb); cmaps["z_bed5"]...)
contour!(axs[1, 2], view(z_bed_t0, i1_wais:i2_wais, j1_wais:j2_wais), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)
contour!(axs[2, 2], view(z_bed_t1, i1_wais:i2_wais, j1_wais:j2_wais), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)
contour!(axs[1, 3], view(z_bed_t0, i1_asb:i2_asb, j1_asb:j2_asb), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)
contour!(axs[2, 3], view(z_bed_t1, i1_asb:i2_asb, j1_asb:j2_asb), levels = [58],
    colormap = cgrad([:gray20, :gray20]), linewidth = lw1)

for (ax, c) in zip(vec([axs[i, j] for j in 2:ncols, i in 1:nrows]), framecolors)
    ax.bottomspinecolor = c
    ax.topspinecolor = c
    ax.leftspinecolor = c
    ax.rightspinecolor = c
    ax.spinewidth = lw2
end

Colorbar(fig[nrows+1, 1:2], label = "Bed elevation (km)", width = Relative(0.5), flipaxis = false,
    vertical = false, ticks = latexifyticks(-4:2, 1e3); cmaps["z_bed5"]...)
elem_1 = LineElement(color = :gray20, linewidth = 2)
Legend(fig[nrows+1, 2:3], [elem_1], [L"$z_b = 58 \, \mathrm{m}$ isoline"], valign = :top,
    framevisible = false)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)
save(plotsdir("16km/bedrock/bed.png"), fig)