include("../../../intro.jl")

hr = HeatmapRtip{Float32}(visc_cases)
dirs = [datadir("output/ais/v2/ramps/wsb/$i") for i in 1:9]
for dir in dirs
    aggregate_xp!(hr, dir)
end
@show [length(f) for f in hr.f]

# Compute values related to 2D maps
k = 3
rsl_index = 7
file_rsl_ref = joinpath(hr.paths[k][rsl_index], "yelmo2Dsm.nc")
x = ncread(file_rsl_ref, "xc")
y = ncread(file_rsl_ref, "yc")
X = ncread(file_rsl_ref, "x2D")
Y = ncread(file_rsl_ref, "y2D")

file_fi_ref = joinpath(hr.paths[k][rsl_index], "fastisostasy.nc")
# x_fi = ncread(file_fi_ref, "xc")
# y_fi = ncread(file_fi_ref, "yc")

nx, ny = size(X)
i1, i2, i3, i4 = [200, 230, 280, 300]
j1, j2, j3, j4 = [40, 55, 70, 140]
x1, x2, x3, x4 = x[i1], x[i2], x[i3], x[i4]
y1, y2, y3, y4 = y[j1], y[j2], y[j3], y[j4]
# ii1, ii2, ii3, ii4 = [findfirst(x_fi .>= x[i]) for i in (i1, i2, i3, i4)]
# jj1, jj2, jj3, jj4 = [findfirst(y_fi .>= y[j]) for j in (j1, j2, j3, j4)]

f_bif = fill(5.85, hr.n_visc_cases)
crange = (39, 45)
trange = (2, 7)
vrange = (47, 52)
zrange = (-430, -410)
rtime_hm_jet = get_rtime_idx(unique(hr.rtime[k]))
rtime_k = get_rtime_idx(hr.rtime[k])

set_theme!(theme_latexfonts())
fs = 20
ms1 = 5
ms2 = 20
axasp = 1.15
lw = 3
rw_cbar = 0.6
cmap = (colormap = cgrad([:lightcoral, :white, :cornflowerblue]), colorrange = crange,
    lowclip = :lightcoral, highclip = :cornflowerblue)
viscmarker = :star5

fig_wsb = Figure(size = (1100, 750), fontsize = fs)
ha = fig_wsb[1, 1:3] = GridLayout()
ra = fig_wsb[2, 1] = GridLayout()
ba = fig_wsb[2, 2:3] = GridLayout()

ax_bed = Axis(ra[1, 1], valign = :top)
ax_vol = Axis(ra[2, 1], valign = :top)
ax_hm = [Axis(ba[1, j], aspect = DataAspect()) for j in 1:2]
axs = [Axis(ha[1, j], aspect = AxisAspect(1)) for j in 1:hr.n_visc_cases]

ax_bed.ylabel = "Mean bed elevation (m)"
xlims!(ax_bed, trange)
ylims!(ax_bed, zrange)
ax_bed.yticks = zrange[1]:10:zrange[2]
ax_bed.yaxisposition = :left
ax_bed.xticklabelsvisible = false
ax_bed.xticksvisible = false

ax_vol.ylabel = L"AIS volume $\mathrm{(m \, SLE)}$    "
ax_vol.xticklabelsvisible = true
ax_vol.xticksvisible = true
xlims!(ax_vol, trange)
ylims!(ax_vol, vrange)
ax_vol.yticks = vrange[1]:1:vrange[2]-1
ax_vol.yaxisposition = :left
ax_vol.xlabel = "Time (kyr)"

axs[1].ylabel = "Max. GMT warming (K)"
axs[hr.n_visc_cases].yaxisposition = :right

for j in eachindex(hr.visc_cases)
    
    if length(hr.f[j]) > 1
        hrtime_idx = get_rtime_idx(hr.rtime[j])
        heatmap!(axs[j], hrtime_idx, hr.f[j] ./ pa .+ f_pd, hr.V[j]; cmap...)
        scatter!(axs[j], hrtime_idx, hr.f[j] ./ pa .+ f_pd, markersize = ms1,
            color = :white)
        scatter!(axs[j], hrtime_idx[rsl_index], hr.f[j][rsl_index] ./ pa .+ f_pd,
            markersize = ms2, color = viscmap[j], marker = viscmarker)
        hlines!(axs[j], [f_bif[j]], color = :gray10, linestyle = :dash, linewidth = 4)
    end
    # vlines!(axs[j], log10.([dfdt_min, dfdt_max]), color = :gray20, linewidth = 3)

    axs[j].title = visc_num_labels[j]
    axs[j].xticks = (rtime_hm_jet,
        reverse(["8e-4", "1e-3", "2e-3", "3e-3", "6e-3", "6e-2", "6e-1", "6e0", "6e1"]))
    vlines!(axs[j], [6.8, 7.8], color = :gray50, linewidth = 3, linestyle = :dash)
    axs[j].xticklabelrotation = π / 2
    ylims!(axs[j], extrema(hr.f[1]) ./ pa .+ f_pd .+ (-0.05, 0.05))

    if 1 < j # < hr.n_visc_cases
        axs[j].yticksvisible = false
        axs[j].yticklabelsvisible = false
    end

    file1D = joinpath(hr.paths[j][rsl_index], "yelmo1D.nc")
    lines!(ax_vol, ncread(file1D, "time") ./ 1e3, ncread(file1D, "V_sle"),
        color = viscmap[j], linewidth = 2, label = visc_labels[j])

    file = joinpath(hr.paths[j][rsl_index], "yelmo2Dsm.nc")
    rsl = vec(mean(
        ncread(file, "z_bed", start = [i2, j2, 1], count = [i3-i2, j3-j2, -1]),
        dims = (1, 2),
    ))
    t2D = ncread(file, "time")
    lines!(ax_bed, t2D ./ 1e3, rsl, color = viscmap[j], linewidth = 2, label = visc_labels[j])
end
axs[3].xlabel = L"GMT warming rate $\mathrm{(K \, century^{-1})}$"
axislegend(ax_vol, position = :rt, nbanks = 3, labelsize = 14)

# The inset axis
dcrop = 10
inset_axs = [
    Axis(fig_wsb[2, i+1],
    aspect = DataAspect(),
    width=Relative(0.3),
    height=Relative(0.3),
    halign=0.02,
    valign=0.97,
    backgroundcolor=:white) for i in [1, 2]]
[hidedecorations!(inset_ax) for inset_ax in inset_axs]

hm_idx = [1, 5]
for i in eachindex(hm_idx)
    
    j = hm_idx[i]
    hidedecorations!(ax_hm[i])
    file = joinpath(hr.paths[j][rsl_index], "yelmo2Dsm.nc")
    time = ncread(file, "time")
    nt = length(time)
    H_ice_ref = ncread(file, "H_ice", start = [1, 1, 1], count = [-1, -1, 1])[:, :, 1]
    H_ice = ncread(file, "H_ice", start = [1, 1, nt], count = [-1, -1, 1])[:, :, 1]
    z_bed = ncread(file, "z_bed", start = [1, 1, nt], count = [-1, -1, 1])[:, :, 1]
    # f_grnd_ref = ncread(file, "f_grnd", start = [i1, i2, 1], count = [d1, d2, 1])[:, :, 1]
    # f_grnd = ncread(file, "f_grnd", start = [i1, i2, nt], count = [d1, d2, 1])[:, :, 1]
    # f_grnd_ref_glob = ncread(file, "f_grnd", start = [1, 1, 1], count = [-1, -1, 1])[:, :, 1]
    heatmap!(ax_hm[i], x, y, z_bed; cmaps["z_bed"]...)
    heatmap!(ax_hm[i], x, y, H_ice + z_bed; cmaps["z_srf"]...)
    xlims!(ax_hm[i], (x1, x4))
    ylims!(ax_hm[i], (y1, y4))

    file1D = joinpath(hr.paths[j][rsl_index], "yelmo1D.nc")
    # contour!(ax_hm[i], f_grnd_ref, color = :black, linewidth = 2, levels = [0.5])
    # contour!(ax_hm[i], f_grnd, color = :red, linewidth = 2, levels = [0.5])
    contour!(ax_hm[i], x, y,
        ((x2 .<= X .<= x3) .&&
        (y2 .<= Y .<= y3)),
        color = :orange,
        linewidth = 2, levels = [0.5])
    heatmap!(inset_axs[i], x, y, H_ice_ref .> 1e-8, colorrange = (1e-8, 1),
        colormap = cgrad([:white, :gray70]), lowclip = :white, highclip = :gray70)
    # contour!(inset_axs[i], f_grnd_ref_glob, color = :black, linewidth = 2, levels = [0.5])
    contour!(inset_axs[i], x, y, ((x1 .<= X .<= x4) .&& (y1 .<= Y .<= y4)),
        color = :black, linewidth = 2)
    xlims!(inset_axs[i], x[dcrop], x[end-dcrop])
    ylims!(inset_axs[i], y[dcrop], y[end-dcrop])
    # text!(inset_axs[i], 10, 5, font = :bold, color = :black, fontsize = fs-4,
    #     text = "t=0 kyr")
        #text = "t = $(Int(round(time[nt] / 1e3, digits = 0))) kyr")
    ax_hm[i].leftspinecolor = viscmap[j]
    ax_hm[i].rightspinecolor = viscmap[j]
    ax_hm[i].topspinecolor = viscmap[j]
    ax_hm[i].bottomspinecolor = viscmap[j]
    ax_hm[i].spinewidth = 5
end

text!(ax_bed, trange[1] + 0.2, zrange[2] - 5, text = "f", color = :black, font = :bold)
text!(ax_vol, trange[1] + 0.2, vrange[1] + 0.3, text = "g", color = :black, font = :bold)
text!(ax_hm[1], x1 + 80, y1 + 60, text = "h", color = :white, font = :bold)
text!(ax_hm[2], x1 + 80, y1 + 60, text = "i", color = :white, font = :bold)

Colorbar(fig_wsb[3, 2], label = L"Final AIS volume $\mathrm{(m \, SLE)}$", vertical = false,
    width = Relative(rw_cbar), flipaxis = false, halign = :left; cmap...)
Colorbar(fig_wsb[3, 2:3], label = "Ice surface elevation (km)", vertical = false,
    width = Relative(rw_cbar/2), flipaxis = false, ticks = (vcat(1, 1e3:1e3:4e3),
    string.(0:4)); cmaps["z_srf"]...)
Colorbar(fig_wsb[3, 3], label = "Bed elevation (km)", vertical = false,
    width = Relative(rw_cbar), flipaxis = false, ticks = latexifyticks(-6:2:4, 1e3),
    halign = :right; cmaps["z_bed"]...)

colsize!(ra, 1, 280)
rowsize!(fig_wsb.layout, 2, 350)
rowgap!(fig_wsb.layout, 10)
rowgap!(ra, 5)
rowgap!(fig_wsb.layout, 2, -40)
colgap!(ha, 5)
fig_wsb

save(plotsdir("v2/rtip/ramp-heatmap-wsb.png"), fig_wsb)
save(plotsdir("v2/rtip/ramp-heatmap-wsb.pdf"), fig_wsb)