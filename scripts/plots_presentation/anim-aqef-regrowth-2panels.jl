include("../intro.jl")

dir = datadir("output/ais/hyster/16km/regrowth/aqef/0")
file1D = "$dir/yelmo1D.nc"
file2D = "$dir/yelmo2Dsm.nc"
file_bsl = "$dir/bsl.nc"

nx, ny = size(ncread(file2D, "x2D"))
X2D, Y2D = ncread(file2D, "x2D"), ncread(file2D, "y2D")
tol = 1e-8
f_to = 0.25

time_1D = ncread(file1D, "time")
dt_1D = time_1D[2] - time_1D[1]
n_1D = length(time_1D)
t_end = maximum(time_1D) - 1e2  # keep some margin because writing out asynchronous

time_2D = ncread(file2D, "time")
dt_2D = time_2D[2] - time_2D[1]
n_2D = length(time_2D)
dt_scale = Int(dt_2D / dt_1D)

hyst_f_now = ncread(file1D, "hyst_f_now")
bsl = ncread(file_bsl, "bsl")
bmb = ncread(file1D, "bmb")
smb = ncread(file1D, "smb")

H_ice = ncread(file2D, "H_ice")
z_bed = ncread(file2D, "z_bed")
uxy_s = ncread(file2D, "uxy_s")
uxy_s[uxy_s .== 0] .= 1e20
z_srf = H_ice .+ z_bed

z_bed_missing = Array{Union{Missing, Float32}}(undef, size(z_bed))
z_bed_missing .= z_bed
# z_bed_missing[R .> 100, :] .= missing

bmb2D = vec(mean(ncread(file2D, "bmb") .* H_ice .> 0, dims = (1,2)))
smb2D = vec(mean(ncread(file2D, "smb") .* H_ice .> 0, dims = (1,2)))
mb2D = vec(mean(ncread(file2D, "mb_net") .* H_ice .> 0, dims = (1,2)))

k = Observable(1)
H_ice_obs = @lift(H_ice[:, :, $k])
z_bed_obs = @lift(z_bed[:, :, $k])
z_srf_obs = @lift(z_srf[:, :, $k])
uxy_s_obs = @lift(log10.(uxy_s[:, :, $k] .+ tol))
k_1D = @lift(dt_scale*($k-1)+1)
time_1D_obs = @lift(time_1D[1:$k_1D] ./ 1e3)
hyst_f_now_obs = @lift(hyst_f_now[1:$k_1D] .* f_to)
bsl_obs = @lift(bsl[1:$k_1D])

set_theme!(theme_latexfonts())
fs = 22
lw = 3
rw = 0.8
tol = 1e-8
tticks_vals_1D = collect(range(0, stop = n_1D, step = 50e3/dt_1D))
tticks_1D = (tticks_vals_1D, latexify.(Int.(tticks_vals_1D .* dt_1D ./ 1e3)))
tticks_vals_2D = collect(range(0, stop = n_2D, step = 20e3/dt_2D))
tticks_2D = (tticks_vals_2D, latexify.(Int.(tticks_vals_2D .* dt_2D ./ 1e3)))
fig = Figure(size = (1500, 800), fontsize = fs)

ncols = 2
axs = [Axis(fig[1, j], aspect = AxisAspect(1)) for j in 1:ncols]
colors = cgrad(:seaborn_colorblind6, range(0, stop = 1, length=6), categorical = true)

axs[1].yticklabelcolor = colors[1]
axs[1].ylabelcolor = colors[1]
axs[1].xlabel = "Time (kyr)"
axs[1].ylabel = L"$\Delta T_\mathrm{ocn}$ (K)"
xlims!(axs[1], (0, n_1D))
ylims!(axs[1], (-0.1, 6))
axs[1].xticks = tticks_1D

ax2 = Axis(fig[1, 1], aspect = AxisAspect(1))
ax2.yaxisposition = :right
ax2.ylabel = "Barystatic sea level contribution (m)"
ax2.yticklabelcolor = colors[2]
ax2.ylabelcolor = colors[2]
hidexdecorations!(ax2)
xlims!(ax2, (0, n_1D))
ylims!(ax2, (-1, 60))

hidedecorations!(axs[ncols])
# hidespines!(axs[ncols])

lines!(axs[1], hyst_f_now_obs, color = colors[1], linewidth = lw)
lines!(ax2, bsl_obs, color = colors[2], linewidth = lw)

zmap = cmaps["z_bed"]
heatmap!(axs[ncols], z_bed_obs; zmap...)
heatmap!(axs[ncols], uxy_s_obs; cmaps["uxy_s3"]...)
contour!(axs[ncols], z_srf_obs, levels = [-1e3, 0, 1e3, 2e3, 3e3, 4e3],
    linewidth = 1, color = :black)

Colorbar(fig[0, ncols], vertical = false, flipaxis = true, width = Relative(rw),
    label = "Bed elevation (km)", ticks = latexifyticks(-4:2, 1e3); zmap...)
Colorbar(fig[2, ncols], vertical = false, flipaxis = false, width = Relative(rw),
    label = "Log10 surface velocity (m/yr)", ticks = latexifyticks(1:4); cmaps["uxy_s3"]...)

rowsize!(fig.layout, 0, 10)
rowsize!(fig.layout, 2, 10)
rowgap!(fig.layout, 1, 10)
rowgap!(fig.layout, 2, -45)
fig

record(fig, plotsdir("16km/anim-regrowth-aqef.mp4"), 1:n_2D, framerate = 24) do i
    k[] = i
end