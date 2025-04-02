include("../intro.jl")

T = Float32
visc_type = "normvisc"    # "equil" or "aqef"
regrowth_dir = datadir("output/ais/hyster/16km/regrowth")
xps = [
    "$regrowth_dir/aqef/pmpt-$visc_type-dpr",
    "$regrowth_dir/aqef/pmpt-$visc_type-upl",
    "$regrowth_dir/aqef/pmpt-$visc_type-fastnormforcing-restarted",
    "$regrowth_dir/aqef/pmpt-$visc_type-fastnormforcing",
]
xp_labels = [
    "DPR",
    "UPL",
    nothing,
    "REF",
]
lws = [3, 3, 5, 5]
cycling_colors = [2, 3, 1, 1]
aqef = AQEFResults(T, xps)

eqldir = datadir("output/ais/hyster/16km/regrowth/equil")
eql = EquilResults(T, eqldir)

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
heatmap_frames = "aqef"    # "equil" or "aqef"
xp_idx = aqef.n_xps
f_ref = aqef.f[end] ./ polar_amplification
stiching_forcing = 2
stiching_idx = findfirst(f_ref .<= stiching_forcing)

set_theme!(theme_latexfonts())
ms1, ms2 = 8, 15
nrows, ncols = 3, 4
forcing_frames = reshape([nothing, 7.6, 7.4, 7, 6.8, 3, 2.2, 2, 1.6, 1.2, 0.2, 0.1], ncols, nrows)'
state_labels = latexify.(reshape(0:11, ncols, nrows)')
fig = Figure(size=(1700, 1320), fontsize = 24)
axs = [Axis(fig[i+1, j], aspect = AxisAspect(1)) for i in 1:nrows, j in 1:ncols]
s = 50
for k in 1:aqef.n_xps
    if k < aqef.n_xps
        lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification, aqef.V_sle[k][1:s:end],
            linewidth = lws[k], label = xp_labels[k], color = Cycled(cycling_colors[k]))
    else
        lines!(axs[1, 1], aqef.f[k][1:s:stiching_idx] ./ polar_amplification,
            aqef.V_sle[k][1:s:stiching_idx], linewidth = lws[k], label = xp_labels[k],
            color = Cycled(cycling_colors[k]))
    end
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
xlims!(axs[1, 1], -2, 11)
fig

file2D_restarted = joinpath(aqef.xps[xp_idx - 1], "0", "yelmo2D.nc")
file2D = joinpath(aqef.xps[xp_idx], "0", "yelmo2D.nc")
X = ncread(file2D, "x2D")
Y = ncread(file2D, "y2D")
xc = ncread(file2D, "xc")
yc = ncread(file2D, "yc")
nx, ny = size(ncread(file2D, "x2D"))

crop = 20
ii = crop+1:nx-crop
jj = crop+1:ny-crop
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
for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    if i > 1 || j > 1
        forcing = forcing_frames[i, j]
        if heatmap_frames == "aqef" || (i == 1 && j == 2)
            
            if forcing > stiching_forcing
                file2D_plot = file2D
                plot_index = xp_idx
            else
                file2D_plot = file2D_restarted
                plot_index = xp_idx - 1
            end

            i3 = argmin( ( aqef.f[plot_index] ./ polar_amplification .-
                forcing) .^ 2)
            f_eq = aqef.f[plot_index][i3] ./ polar_amplification
            V_eq = aqef.V_sle[plot_index][i3]
            frame_index = argmin(abs.(aqef.t_2D[plot_index] .- aqef.t_1D[plot_index][i3]))
            z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D_plot,
                var_names_2D, frame_index)
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

        if mod(j, 2) == 1
            if i == 1 && j == 3
                vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.2, color = :gray, label = "Bifurcations")
            else
                vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.2, color = :gray)
            end
        end

        scatter!(axs[1, 1], f_eq, V_eq, color = :red, markersize = ms2)
        text!(axs[1, 1], f_eq, V_eq, text = state_labels[i, j],
            color = :grey10, fontsize = 30, font = :bold, offset = text_offsets[i, j])

        hidedecorations!(axs[i, j])
        heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
        heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
        contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
            color = :red, linewidth = 2)
        text!(axs[i, j], -2500, -2500, color = :white, font = :bold,
            text=state_labels[i, j], fontsize = 30)
        xlims!(axs[i, j], extrema(XX))
        ylims!(axs[i, j], extrema(YY))
    end
end

axislegend(axs[1, 1], position = :rt, nbanks = 1)
relwidth = 0.8
Colorbar(fig[1, 2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[1, 3], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
Legend(fig[1, 4], [elem_1], ["Grounding line"])

rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -20)
rowsize!(fig.layout, 1, 1)
rowsize!(fig.layout, 2, 400)
rowsize!(fig.layout, 3, 400)
rowsize!(fig.layout, 4, 400)
colsize!(fig.layout, 1, 400)
colsize!(fig.layout, 2, 400)
colsize!(fig.layout, 3, 400)
colsize!(fig.layout, 4, 400)
save(plotsdir("16km/hysteresis/regrowth.png"), fig)