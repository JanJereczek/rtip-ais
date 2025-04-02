include("../intro.jl")

file = datadir("output/ais/hyster/16km/retreat/aqef/pmpt-normvisc-normforcing-withrestarts/0/yelmo2D.nc")
t2D = ncread(file, "time")
z_bed_ref = ncslice(file, "z_bed", 1)
z_bed_end = ncslice(file, "z_bed", length(t2D))

j2020 = "/p/projects/megarun/ice_data/ISMIP6/Antarctica/ANT-16KM/Ocean/ANT-16KM_OCEAN_ISMIP6_J20.nc"
x, y, z = ncread(j2020, "x"), ncread(j2020, "y"), ncread(j2020, "z")
z_sorted = sort(z)
tf = ncread(j2020, "tf")
tf = reverse(tf, dims = 3)
itp = linear_interpolation((x, y, z_sorted), tf, extrapolation_bc = NaN)
X, Y = ndgrid(x, y)
tf_ref = itp.(X, Y, z_bed_ref)
tf_end = itp.(X, Y, z_bed_end)

set_theme!(theme_latexfonts())
nrows, ncols = 2, 3
fig = Figure(size = (1200, 900), fontsize = 24)
axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(axs[i, j]) for i in 1:nrows, j in 1:ncols]

heatmap!(axs[1, 1], z_bed_ref; cmaps["z_bed"]...)
heatmap!(axs[1, 2], z_bed_end; cmaps["z_bed"]...)
heatmap!(axs[1, 3], z_bed_end .- z_bed_ref; cmaps["dz_bed"]...)

heatmap!(axs[2, 1], tf_ref; cmaps["tf"]...)
heatmap!(axs[2, 2], tf_end; cmaps["tf"]...)
heatmap!(axs[2, 3], tf_end .- tf_ref; cmaps["dtf"]...)

Colorbar(fig[nrows + 1, 1], vertical = false, flipaxis = false, width = Relative(0.8),
    label = "Bed elevation (km)", ticks = latexifyticks(-4:2:2, 1e3); cmaps["z_bed"]...)
Colorbar(fig[nrows + 1, 2], vertical = false, flipaxis = false, width = Relative(0.8),
    label = "Thermal forcing (K)"; cmaps["tf"]...)
Colorbar(fig[nrows + 1, 3], vertical = false, flipaxis = false, width = Relative(0.8),
    label = "Thermal forcing difference (K)"; cmaps["dtf"]...)

axs[1, 1].title = "Before uplift"
axs[1, 2].title = "After uplift"
axs[1, 3].title = "After minus Before"

save(plotsdir("16km/bed-thermal-uplifted.png"), fig)