include("../intro.jl")

T = Float32
visc_type = "lowvisc"    # "equil" or "aqef"
heatmap_frames = "equil"    # "equil" or "aqef"

xps = [
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-normforcing-withrestarts"),
]
xp_labels = [
    "REF",
]
lws = [5]
aqef = AQEFResults(T, xps)
eqldir = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
eql = EquilResults(T, eqldir)

############################################################################
# A simple mb plot for the appendix
############################################################################
plot_mb = false

if plot_mb
    lw = 3
    set_theme!(theme_latexfonts())
    fig = Figure(size = (800, 500), fontsize = 24)

    ax2 = Axis(fig[1, 1], yaxisposition = :right)
    hidexdecorations!(ax2)
    dvdt = diff(aqef.V_sle[1]) ./ diff(aqef.t_1D[1])
    s = 10
    lines!(ax2, aqef.t_1D[1][2:s:end] ./ 1f3, dvdt[1:s:end] .* 1f3, label = "Volume change",
        linewidth = lw, color = :gray70)
    axislegend(ax2, position = :rb, labelsize = 20)
    xlims!(ax2, 0, 200)
    ylims!(ax2, -2, 2)
    ax2.ylabel = "Volume change (mmSLE/yr)"
    ax2.yticks = -2:1:2

    ax = Axis(fig[1, 1])
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.bmb[1], label = "basal", linewidth = lw)
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.smb[1], label = "surface", linewidth = lw)
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.cmb[1], label = "calving", linewidth = lw)
    lines!(ax, aqef.t_2D[1] ./ 1f3, aqef.bmb[1] .+ aqef.smb[1] .+ aqef.cmb[1], 
    linewidth = lw, label = "total")
    axislegend(ax, position = :rt, labelsize = 20)
    xlims!(ax, 0, 200)
    ylims!(ax, -0.1, 0.1)
    ax.xlabel = "Time (kyr)"
    ax.ylabel = "Mass balance (m/yr)"
    ax.yticks = -0.1:0.05:0.1

    save(plotsdir("16km/hysteresis/fluxes.png"), fig)
end

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
xp_idx = aqef.n_xps

cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)
set_theme!(theme_latexfonts())
ms1, ms2 = 8, 15
nrows, ncols = 3, 4
forcing_frames = reshape([0, 0, 1.15, 1.35, 4.2, 4.4, 5.8, 6.0, 6.7, 7.1, 7.6, 7.8], ncols, nrows)'
state_labels = latexify.(reshape(0:11, ncols, nrows)')
fig = Figure(size=(1400, 1050), fontsize = 24)
axs = [Axis(fig[i+1, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]
s = 50
for k in 1:aqef.n_xps
    lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k])
end
scatterlines!(axs[1, 1], eql.f ./ polar_amplification, eql.V_sle;
    linewidth = lws[xp_idx], color = :black, label = "EQL", markersize = ms1)
axs[1, 1].xticks = 0:2:12
axs[1, 1].xminorticks = 0:0.2:12
axs[1, 1].yticks = 0:10:60
axs[1, 1].yminorticks = IntervalsBetween(10)
axs[1, 1].xminorgridvisible = true
axs[1, 1].yminorgridvisible = true
axs[1, 1].xaxisposition = :top
axs[1, 1].xlabel = L"GMT anomaly $f$ (K)"
axs[1, 1].ylabel = L"AIS volume $V_\mathrm{af}$ (m SLE)"
ylims!(axs[1, 1], 0, 60)
xlims!(axs[1, 1], 0, 11)
fig

# equil_opts = (label = "Equilibrium", color = :black, linewidth = 4)
# scatterlines!(axs[1, 1], f_equil, V_equil; equil_opts...)

file2D = joinpath(aqef.xps[xp_idx], "0", "yelmo2D.nc")
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))

ii = cropx+1:nx-cropx
jj = cropy+1:ny-cropy
XX = X[ii, jj]
YY = Y[ii, jj]

text_offsets = permutedims(reshape([
    (0, 0),         # 0
    (7, -33),         # 1
    (7, -12),      # 2
    (-20, -32),     # 3
    (7, -10),       # 4
    (-17, -10),       # 5
    (7, -11),     # 6
    (-18, -20),       # 7
    (7, -15),       # 8
    (-22, -15),     # 9
    (5, -10),       # 10
    (10, -22),       # 11
], ncols, nrows))

var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]

region_labels = ["WAIS", "EAIS"]
region_positions = [(-1300, -300), (200, 0)]
subregion_labels = ["ASE", "WSB", "RSB", "ASB"]
subregion_positions = [(-1450, -600), (600, -1750), (-400, 900), (1600, -800)]
region_color = :black
subregion_color = :honeydew
region_fontsize = 28
subregion_fontsize = 22

for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    if i > 1 || j > 1
        forcing = forcing_frames[i, j]

        if mod(j, 2) == 1
            if i == 1 && j == 3
                vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.2,
                    color = :gray, label = "Bifurcation")
            else
                vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.2,
                    color = :gray)
            end
        end

        if heatmap_frames == "aqef" || (i == 1 && j == 2)
            i3 = findfirst(aqef.f[xp_idx] ./ polar_amplification .>= forcing)
            f_eq, V_eq = aqef.f[xp_idx][i3] ./ polar_amplification, aqef.V_sle[xp_idx][i3]
            frame_index = argmin(abs.(aqef.t_2D[xp_idx] .- aqef.t_1D[xp_idx][i3]))
            z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D,
                frame_index)
        else
            i3 = argmin((eql.f ./ polar_amplification .- forcing) .^ 2)
            file_2D_equil = datadir("$eqldir/$(string(i3))/0/yelmo2D.nc")
            @show file_2D_equil
            f_eq, V_eq = eql.f[i3] ./polar_amplification, eql.V_sle[i3]
            frame_index = length(ncread(file_2D_equil, "time"))
            z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file_2D_equil,
                var_names_2D, frame_index)
        end
        @show i3, f_eq, V_eq

        scatter!(axs[1, 1], f_eq, V_eq, color = :red, markersize = ms2)
        text!(axs[1, 1], f_eq, V_eq, text = state_labels[i, j],
            color = :grey10, fontsize = 30, font = :bold, offset = text_offsets[i, j])

        hidedecorations!(axs[i, j])
        heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
        heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
        contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
            color = :red, linewidth = 2)
        # if transects[i, j] !== nothing
        #     lines!(axs[i, j], transects[i, j].x, transects[i, j].y; color = :orange, linewidth = 3)
        # end
        text!(axs[i, j], -2500, -2500, color = :white, font = :bold, text=state_labels[i, j],
            fontsize = 30)
        xlims!(axs[i, j], extrema(XX))
        ylims!(axs[i, j], extrema(YY))
    end
end

for i in eachindex(region_labels)
    text!(axs[1, 2], region_positions[i]..., text = region_labels[i], color = region_color,
        fontsize = region_fontsize, font = :bold)
end
for i in eachindex(subregion_labels)
    text!(axs[1, 2], subregion_positions[i]..., text = subregion_labels[i], color = subregion_color,
        fontsize = subregion_fontsize, font = :bold)
end

axislegend(axs[1, 1], position = :lb, nbanks = 1)
relwidth = 0.8
Colorbar(fig[1, 2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[1, 3], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :orange, linewidth = 3)
Legend(fig[1, 4], [elem_1], ["Grounding line"])

rowsize_base = 300
rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -20)
rowsize!(fig.layout, 1, 1)
rowsize!(fig.layout, 2, rowsize_base)
rowsize!(fig.layout, 3, rowsize_base)
rowsize!(fig.layout, 4, rowsize_base)
colsize!(fig.layout, 1, rowsize_base*aratio)
colsize!(fig.layout, 2, rowsize_base*aratio)
colsize!(fig.layout, 3, rowsize_base*aratio)
colsize!(fig.layout, 4, rowsize_base*aratio)
save(plotsdir("16km/hysteresis/min-retreat-$visc_type.png"), fig)
save(plotsdir("16km/hysteresis/min-retreat-$visc_type.pdf"), fig)