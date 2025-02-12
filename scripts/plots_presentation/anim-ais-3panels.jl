include("../intro.jl")

file1D = datadir("output/ais/hyster/aqef/20K/1/yelmo1D.nc")
file2D = datadir("output/ais/hyster/aqef/20K/1/yelmo2Dsm.nc")

nx, ny = 191, 191
x = 1:nx
y = 1:ny
X = repeat(x, 1, ny)
Y = repeat(y', nx, 1)
R = sqrt.((X .- nx/2).^2 + (Y .- ny/2).^2)

f_to = 0.25
t_end = 1.9e5
dt_1D = 5
dt_2D = 200
dt_scale = Int(dt_2D / dt_1D)
time_1D = collect(0:dt_1D:t_end)
n_1D = length(time_1D)
n_1D_wrong = length(ncread(file1D, "time"))

stride_1D = 2000
time_1D = ncread(file1D, "time")[1:stride_1D:end]
hyst_f_now = ncread(file1D, "hyst_f_now")[1:stride_1D:end]
bsl = ncread(file1D, "bsl")[1:stride_1D:end]
bmb = ncread(file1D, "bmb")[1:stride_1D:end]
smb = ncread(file1D, "smb")[1:stride_1D:end]

time_2D = ncread(file2D, "time")
n_2D = length(time_2D)
H_ice = ncread(file2D, "H_ice")
z_bed = ncread(file2D, "z_bed")
uxy_s = ncread(file2D, "uxy_s")
uxy_s[uxy_s .== 0] .= 1e20
z_srf = H_ice .+ z_bed

z_bed_missing = Array{Union{Missing, Float32}}(undef, size(z_bed))
z_bed_missing .= z_bed
z_bed_missing[R .> 100, :] .= missing

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
bmb_obs = @lift(bmb[1:$k_1D])
smb_obs = @lift(smb[1:$k_1D])
bmb2D_obs = @lift(bmb2D[1:$k])
smb2D_obs = @lift(smb2D[1:$k])
mb2D_obs = @lift(cmb2D[1:$k])

set_theme!(theme_latexfonts())
fs = 22
lw = 3
rw = 0.8
tol = 1e-8
tticks_vals_1D = collect(range(0, stop = n_1D, step = 20e3/dt_1D))
tticks_1D = (tticks_vals_1D, latexify.(Int.(tticks_vals_1D .* dt_1D ./ 1e3)))
tticks_vals_2D = collect(range(0, stop = n_2D, step = 20e3/dt_2D))
tticks_2D = (tticks_vals_2D, latexify.(Int.(tticks_vals_2D .* dt_2D ./ 1e3)))
fig = Figure(size = (1500, 600), fontsize = fs)

axs = [Axis(fig[1, j], aspect = AxisAspect(1)) for j in 1:3]
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

axs[2].xlabel = "Time (kyr)"
axs[2].ylabel = "Mass balance (m/yr)"
axs[2].yaxisposition = :right

hidedecorations!(axs[3])
# hidespines!(axs[3])

lines!(axs[1], hyst_f_now_obs, color = colors[1], linewidth = lw)
lines!(ax2, bsl_obs, color = colors[2], linewidth = lw)

mb = "1D"
if mb == "1D"
    lines!(axs[2], smb_obs, color = colors[3], linewidth = lw, label = "surface")
    lines!(axs[2], bmb_obs, color = colors[4], linewidth = lw, label = "basal")
    ylims!(axs[2], (-0.2, 0.3))
    xlims!(axs[2], (0, n_1D))
    axs[2].xticks = tticks_1D

elseif mb == "2D"
    lines!(axs[2], smb2D_obs, color = colors[3], linewidth = lw, label = "surface")
    lines!(axs[2], bmb2D_obs, color = colors[4], linewidth = lw, label = "basal")
    lines!(axs[2], mb2D_obs, color = colors[5], linewidth = lw, label = "net")
    xlims!(axs[2], (0, n_2D))
    ylims!(axs[2], (-0.1, 0.1))
    axs[2].xticks = tticks_2D
end
axislegend(axs[2], position = :lc)

heatmap!(axs[3], z_bed_obs; cmaps["z_bed3"]...)
heatmap!(axs[3], uxy_s_obs; cmaps["uxy_s2"]...)
contour!(axs[3], z_srf_obs, levels = [-1e3, 0, 1e3, 2e3, 3e3, 4e3],
    linewidth = 1, color = :black)

Colorbar(fig[0, 3], vertical = false, flipaxis = true, width = Relative(rw),
    label = "Bed elevation (km)", ticks = latexifyticks(-1:2, 1e3); cmaps["z_bed3"]...)
Colorbar(fig[2, 3], vertical = false, flipaxis = false, width = Relative(rw),
    label = "Log10 surface velocity (m/yr)"; cmaps["uxy_s2"]...)

rowsize!(fig.layout, 0, 10)
rowsize!(fig.layout, 2, 10)
rowgap!(fig.layout, 1, 5)
rowgap!(fig.layout, 2, -45)
fig

record(fig, "anim-ais.mp4", 1:n_2D, framerate = 24) do i
    k[] = i
end