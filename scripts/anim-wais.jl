include("intro.jl")

file = datadir("output/ais/hyster/aqef/20K/3/yelmo2Dwais.nc")

vars = ["H_ice", "z_bed", "f_grnd", "uxy_s", "z_sl", "T_shlf", "smb"]
vars = [ncread(file, v) for v in vars]
t = ncread(file, "time")

nrows, ncols = 2, 3
fig = Figure(size=(1600, 1280), fontsize = 30)
axs = [Axis(fig[i, j], aspect = DataAspect()) for i in 1:nrows, j in 1:ncols]
[hidedecorations!(ax) for ax in axs]

nt = size(vars[1], 3)
k = Observable(nt)
heatmap!(axs[1, 1], @lift(vars[2][:, :, $k]); cmaps["z_bed"]...)
heatmap!(axs[1, 1], @lift((vars[1][:, :, $k] .+ vars[2][:, :, $k]) .*
    (vars[1][:, :, $k] .> 1)); cmaps["z_srf"]...)
heatmap!(axs[1, 2], @lift(log10.(vars[4][:, :, $k] .+ 1e-8)); cmaps["uxy_s"]...)
contour!(axs[1, 2], @lift(vars[3][:, :, $k] .* (vars[1][:, :, $k] .> 1)); cmaps["f_grnd"]...)
heatmap!(axs[1, 3], @lift(vars[5][:, :, $k]); cmaps["z_sl"]...)
heatmap!(axs[2, 3], @lift(vars[2][:, :, $k] - vars[2][:, :, 1]); cmaps["u_bed"]...)
heatmap!(axs[2, 1], @lift(vars[6][:, :, $k]); cmaps["T_shlf"]...)
heatmap!(axs[2, 2], @lift(vars[7][:, :, $k]); cmaps["smb"]...)

relwidth1 = 0.49
relwidth2 = 0.7

Colorbar(fig[0, 1], vertical = false, width = Relative(relwidth1),
    label = L"$z_\mathrm{bed}$ (km)", ticks = (-4e3:2e3:2e3, latexify.(-4:2:2)), halign = :left;
    cmaps["z_bed"]...)
Colorbar(fig[0, 1], vertical = false, width = Relative(relwidth1),
    label = L"$z_\mathrm{srf}$ (km)", ticks = (1e1:1e3:4e3, latexify.(1:1:4)), halign = :right;
    cmaps["z_srf"]...)

Colorbar(fig[0, 2], vertical = false, width = Relative(relwidth1),
    label = L"$u_\mathrm{s} \: \mathrm{(m \, yr^{-1})}$ ", halign = :left,
    ticks = (0:4, latexify.(["10^$n" for n in 0:4])); cmaps["uxy_s"]...)
# Legend(fig[0, 2], [LineElement(color = :red, linewidth = 3, linestyle = nothing)], [L"gr. line $\,$"], halign = :right, framewidth = 0.5)

Colorbar(fig[0, 3], vertical = false, width = Relative(relwidth),
    label = L"$z_\mathrm{ss}$ (m)", ticks = (-10:5:10, latexify.(-10:5:10)); cmaps["z_sl"]...)

Colorbar(fig[nrows+1, 1], vertical = false, width = Relative(relwidth), flipaxis = false,
    label = L"$T_\mathrm{shlf}$ (K)", ticks = (260:5:280, latexify.(260:5:280)); cmaps["T_shlf"]...)
Colorbar(fig[nrows+1, 2], vertical = false, width = Relative(relwidth), flipaxis = false,
    label = L"$\mathrm{SMB \:  (m \, yr^{-1}})$", ticks = (-60:20:20, latexify.(-60:20:20)); cmaps["smb"]...)
Colorbar(fig[nrows+1, 3], vertical = false, width = Relative(relwidth), flipaxis = false,
    label = L"$u_\mathrm{bed} \: (\times 100 \, \mathrm{m})$ ", ticks = (vcat(-100, 0:200:800), latexify.(vcat(-1, 0:2:8))); cmaps["u_bed"]...)

text!(axs[1, 2], 2, 2, text = @lift("t = $(t[$k]) yr"), color = :black, fontsize = 24)
fig

record(fig, "anim-wais.mp4", 1:nt, framerate = 24) do i
    k[] = i
end