include("../../intro.jl")

file = datadir("BedMachineAntarctica-v3.nc")
z_bed = reverse(ncread(file, "bed"), dims = 2)
h = reverse(ncread(file, "thickness"), dims = 2)

######################################################
# Friction
######################################################

@kwdef struct PowerLaw{T<:AbstractFloat}
    q::T = 0.75
    u0::T = 100
end

@kwdef struct RegularizedCoulombLaw{T<:AbstractFloat}
    q::T = 0.2
    u0::T = 100
end

friction(u::T, law::PowerLaw{T}) where T = u / ( law.u0^law.q + u^(1-law.q) )
friction(u::T, law::RegularizedCoulombLaw{T}) where T = (u / (law.u0 + u))^law.q

u_vec = collect(range(0, stop = 100, step = 0.1))
T = eltype(u_vec)
pl = PowerLaw{T}()
rcl = RegularizedCoulombLaw{T}()

f_power = map(u -> friction(u, pl), u_vec)
f_coulomb = map(u -> friction(u, rcl), u_vec)

set_theme!(theme_latexfonts())
fig = Figure(size = (1400, 580), fontsize = 22)
axf = Axis(fig[1, 1], aspect = AxisAspect(1))
axf.xlabel = L"Sliding velocity $u_b$ ($\mathrm{m \, yr^{-1}}$)"
axf.ylabel = L"Sliding factor $f \left( u_b \right)$ (1)"
lines!(axf, u_vec, f_power, label = L"$\frac{u}{u_0^q + u^{(1-q)}}$ (Garbe et al., 2020)",
    linewidth = 3)
lines!(axf, u_vec, f_coulomb, label = L"$\left( \frac{u}{u_0 + u} \right)^q$ (present study)", linewidth = 3)
axislegend(axf, position = :lt)
fig

######################################################
# TPMP
######################################################

h_af = h * 910 / 1028 + z_bed
nx, ny = size(z_bed)

di, dj = 60, 50
iii1 = 2000:10:5000
jjj1 = 5000:10:8000
iii2 = 9000:10:12000
jjj2 = 7000:10:10000
X1, Y1 = ndgrid(iii1, jjj1)
X2, Y2 = ndgrid(iii2, jjj2)

# thwaites
ii1 = 3520:3670
jj1 = 5680:5830

# amery
ii2 = 10000:10150
jj2 = 8000:8150

# # ross
# ii3 = 4500:8000
# jj3 = 3500:7000

axs = [Axis(fig[1, j], aspect = DataAspect()) for j in 2:3]
inset_axs = [Axis(fig[1, j], aspect = DataAspect(), width=Relative(0.3), height=Relative(0.3),
    halign=0.0, valign=0.09, backgroundcolor=:white) for j in 2:3]
hidedecorations!.([axs..., inset_axs...])
function plot_grounding_zone!(ax, inset_ax, zb, haf, ii, jj, iii, jjj, X, Y, lw, lc)
    heatmap!(ax, view(zb, ii, jj); cmaps["z_bed6"]...)
    contour!(ax, view(haf, ii, jj); levels = [0, 100], color = [:orange, :red],
        linewidth = 3)

    i1 = di
    j1 = dj
    lines!(ax, [i1, i1+10], [j1, j1]; color = lc, linewidth = lw)
    lines!(ax, [i1, i1], [j1-2, j1+2]; color = lc, linewidth = lw)
    lines!(ax, [i1+10, i1+10], [j1-2, j1+2]; color = lc, linewidth = lw)
    text!(ax, i1-2, j1-10, text = "5 km", color = lc)


    heatmap!(inset_ax, view(h, iii, jjj) .> 1e-8, colormap = cgrad([:gray70, :gray95], [0, 1]),
        colorrange = (0, 1), lowclip = :transparent, highclip = :white)
    contour!(inset_ax, (ii[1] .< X .< ii[end]) .&& (jj[1] .< Y .< jj[end]),
        levels = [0.5], color = :red, linewidth = 3)
end

lw, lc = 4, :black
plot_grounding_zone!(axs[1], inset_axs[1], z_bed, h_af, ii1, jj1, iii1, jjj1, X1, Y1, lw, lc)
plot_grounding_zone!(axs[2], inset_axs[2], z_bed, h_af, ii2, jj2, iii2, jjj2, X2, Y2, lw, lc)

elem_1 = LineElement(color = :orange, linewidth = 2)
elem_2 = LineElement(color = :red, linewidth = 2)
Legend(fig[2, 2], [elem_1, elem_2], ["Downstream limit", "Upstream limit"],
    nbanks = 1, valign = :top)
Colorbar(fig[2, 3], vertical = false, label = "Bed elevation (km)", flipaxis = false,
    width = Relative(0.5), valign = :bottom; cmaps["z_bed6"]...)
axf.title = "(a) Friction laws comparison"
axs[1].title = "(b) Thwaites grounding zone"
axs[2].title = "(c) Amery grounding zone"

base = 400
colsize!(fig.layout, 1, base)
colsize!(fig.layout, 2, base)
colsize!(fig.layout, 3, base)
rowsize!(fig.layout, 2, 0)
rowgap!(fig.layout, 1, -60)
fig
save(plotsdir("16km/hysteresis/figA1.png"), fig)
save(plotsdir("16km/hysteresis/figA1.pdf"), fig)