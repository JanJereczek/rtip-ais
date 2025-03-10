include("../intro.jl")
file = datadir("BedMachineAntarctica-v3.nc")
z_bed = reverse(ncread(file, "bed"), dims = 2)
h = reverse(ncread(file, "thickness"), dims = 2)
h_af = h * 910 / 1028 + z_bed
nx, ny = size(z_bed)

i1, j1 = 150, 120
iii = 2000:10:12000
jjj = 3000:10:11000
region = "ross"

if region == "thwaites"
    ii = 3450:3700
    jj = 5700:6000
elseif region == "ross"
    ii = 4500:8000
    jj = 3500:7000
end

set_theme!(theme_latexfonts())
fig = Figure(size = (1300, 1500), fontsize = 32)
ax = Axis(fig[1, 1], aspect = DataAspect())
hidedecorations!(ax)
heatmap!(ax, view(z_bed, ii, jj); cmaps["z_bed6"]...)
# heatmap!(ax, view((z_bed .+ h) .* Float32.(h .> 1e-8), ii, jj); cmaps["z_srf"]...)
contour!(ax, view(h_af, ii, jj); levels = [0, 100], color = [:orange, :red],
    linewidth = 3)

lw = 4
lc = :black
lines!(ax, [i1, i1+10], [j1, j1]; color = lc, linewidth = lw)
lines!(ax, [i1, i1], [j1-2, j1+2]; color = lc, linewidth = lw)
lines!(ax, [i1+10, i1+10], [j1-2, j1+2]; color = lc, linewidth = lw)
text!(ax, i1-2, j1+2, text = "5 km", color = lc)

X, Y = ndgrid(iii, jjj)
inset_ax = Axis(fig[1, 1], aspect = DataAspect(), width=Relative(0.3), height=Relative(0.3),
    halign=0.05, valign=0.03, backgroundcolor=:white)
hidedecorations!(inset_ax)
heatmap!(inset_ax, view(h, iii, jjj) .> 1e-8, colormap = cgrad([:gray70, :gray95], [0, 1]),
    colorrange = (0, 1), lowclip = :transparent, highclip = :white)
contour!(inset_ax, (ii[1] .< X .< ii[end]) .&& (jj[1] .< Y .< jj[end]),
    levels = [0.5], color = :red, linewidth = 3)


save(plotsdir("tpmp-$region.png"), fig)