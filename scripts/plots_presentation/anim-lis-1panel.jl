include("../intro.jl")

sim = "fast"
file2D = datadir("output/yelmo2D_reduced_$sim.nc")

# nx, ny = size(ncread(file2D, "x2D"))
# X2D, Y2D = ncread(file2D, "x2D"), ncread(file2D, "y2D")

time_2D = ncread(file2D, "time")
dt_2D = time_2D[2] - time_2D[1]
n_2D = length(time_2D)

H_ice = ncread(file2D, "H_ice")
z_bed = ncread(file2D, "z_bed")
uxy_s = ncread(file2D, "uxy_s")
# uxy_s[uxy_s .== 0] .= 1e20
z_srf = H_ice .+ z_bed

z_bed_missing = Array{Union{Missing, Float32}}(undef, size(z_bed))
z_bed_missing .= z_bed
# z_bed_missing[R .> 100, :] .= missing

k = Observable(1)
# k_1D = @lift(dt_scale*($k-1)+1)
H_ice_obs = @lift(H_ice[:, :, $k])
z_bed_obs = @lift(z_bed[:, :, $k])
z_srf_obs = @lift(z_srf[:, :, $k])
# uxy_s_obs = @lift(log10.(uxy_s[:, :, $k] .+ tol))
lin_uxy_s_obs = @lift(uxy_s[:, :, $k])
# time_1D_obs = @lift(time_1D[1:$k_1D] ./ 1e3)
# hyst_f_now_obs = @lift(hyst_f_now[1:$k_1D] .* f_to)
t_obs = @lift("t = $(round(time_2D[$k])) yr")

set_theme!(theme_latexfonts())
fs = 22
lw = 3
rw = 0.45
tol = 1e-8
# tticks_vals_2D = collect(range(0, stop = n_2D, step = 20e3/dt_2D))
# tticks_2D = (tticks_vals_2D, latexify.(Int.(tticks_vals_2D .* dt_2D ./ 1e3)))
fig = Figure(size = (900, 800), fontsize = fs)

ax = Axis(fig[1, 1], aspect = AxisAspect(1))
hidedecorations!(ax)
zmap = cmaps["z_bed"]
uxy_s_map4 = cgrad([:white, :darkred], range(0, stop = 1, length = 11), categorical = true)
lin_umap = (colormap = uxy_s_map4, colorrange = (1f-8, 1000),
lowclip = :transparent, highclip = uxy_s_map4[end])
heatmap!(ax, z_bed_obs; zmap...)
# heatmap!(ax, uxy_s_obs; cmaps["uxy_s3"]...)
heatmap!(ax, lin_uxy_s_obs; lin_umap...)
contour!(ax, z_srf_obs, levels = [-1e3, 0, 1e3, 2e3, 3e3, 4e3],
    linewidth = 1, color = :black)

Colorbar(fig[1, 2], vertical = true, flipaxis = true, height = Relative(rw), valign = :top,
    label = "Bed elevation (km)", ticks = latexifyticks(-4:2, 1e3); zmap...)
# Colorbar(fig[1, 2], vertical = true, flipaxis = true, height = Relative(rw), valign = :bottom,
#     label = "Log10 surface velocity (m/yr)", ticks = latexifyticks(1:4); lin_umap...)
Colorbar(fig[1, 2], vertical = true, flipaxis = true, height = Relative(rw), valign = :bottom,
    label = "Surface velocity (m/yr)", ticks = 0:200:1000; lin_umap...)
colgap!(fig.layout, 5)
text!(ax, 2, 205, text = t_obs, color = :white, fontsize = fs, font = :bold)
fig

record(fig, plotsdir("16km/lis-$sim.mp4"), 1:n_2D, framerate = 10) do i
    k[] = i
end