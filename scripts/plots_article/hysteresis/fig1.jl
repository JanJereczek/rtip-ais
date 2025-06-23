include("../../intro.jl")

T = Float32
visc_type = "lowvisc"    # "equil" or "aqef"
heatmap_frames = "aqef"    # "equil" or "aqef"

xps = [
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-ocnforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-normvisc-fastatmforcing"),
    datadir("output/ais/hyster/16km/retreat/aqef/pmpt-$visc_type-normforcing-withrestarts"),
]
xp_labels = [
    "OCN",
    "ATM",
    "REF",
]
aqef = AQEFResults(T, xps)
eqldir = datadir("output/ais/hyster/16km/retreat/equil/pmpt-normvisc-normforcing")
eql = EquilResults(T, eqldir)
lws = [3, 3, 6]
cycling_colors = [
    xpcolors["OCN"],
    xpcolors["ATM"],
    xpcolors["REF"],
]

############################################################################
# The retreat comparison
############################################################################

polar_amplification = 1.8
f_to = 0.25
xp_idx = aqef.n_xps
f2015 = 1.2

cropx, cropy = 20, 35
aratio = (381 - 2*cropx) / (381 - 2*cropy)
set_theme!(theme_latexfonts())
ms1, ms2 = 8, 18
nrows, ncols = 3, 4
forcing_frames = transpose(reshape([0, 0, 1.15, 1.35, 4.55, 4.75, 6.05, 6.25, 6.95, 7.15,
    7.75, 7.95], ncols, nrows))
state_labels = ["a" "b" "c" "d";
    "e" "f" "g" "h";
    "i" "j" "k" "l"]
fig = Figure(size=(1400, 1050), fontsize = 24)
axs = [Axis(fig[i+1, j], aspect = AxisAspect(aratio)) for i in 1:nrows, j in 1:ncols]
s = 400

shade = [(1.2, 1.3), (4.6, 4.7), (6.1, 6.2), (7.0, 7.1), (7.8, 7.9)]
lightshade = [(4.6, 4.7), (8.5, 8.6), (9.9, 10), (10.2, 10.3)]
for i in eachindex(shade)
    vlines!(axs[1, 1], (shade[i][1]:0.01:shade[i][2]) .+ f2015, alpha = 0.2,
        linewidth = 3, color = :gray70)
end
vlines!(axs[1, 1], 1f6, alpha = 0.9, color = :gray70, linewidth = 6, label = "Bifurcation")

text!(axs[1, 1], 10.2, 50, text = "(a)", color = :grey10, fontsize = 30, font = :bold)
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
xlims!(axs[1, 1], 0, 12)
fig

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

xl = [
    (-1900, -900),
    (700, 1400),
    (-600, 400),
    (400, 1400),
    (1300, 2500)
]
yl = [
    (-800, 200),
    (-2300, -1600),
    (800, 1800),
    (-1600, -600),
    (-1200, 0),
]

xlims_frames = permutedims(reshape([nothing, nothing, xl[1], xl[1], xl[2], xl[2],
    xl[3], xl[3], xl[4], xl[4], xl[5], xl[5]], ncols, nrows))
ylims_frames = permutedims(reshape([nothing, nothing, yl[1], yl[1], yl[2], yl[2],
    yl[3], yl[3], yl[4], yl[4], yl[5], yl[5]], ncols, nrows))

text_offsets = permutedims(reshape([
    (0, 0),         # 1
    (-20, -40),       # 2
    (7, -16),       # 3
    (-15, -45),     # 4
    (10, -20),       # 5
    (-20, -45),     # 6
    (7, -20),       # 7
    (-22, -33),     # 8
    (15, -15),       # 9
    (-25, -20),     # 10
    (10, -15),       # 11
    (-25, -20),     # 12
], ncols, nrows))

var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]

region_labels = ["WAIS", "EAIS"]
region_positions = [(-1300, -300), (200, 0)]
subregion_labels = ["ASE", "WSB", "RSB", "ASB"]
subregion_positions = [(-1450, -600), (600, -1750), (-450, 900), (1600, -800)]
region_color = :black
subregion_color = :white
region_fontsize = 28
subregion_fontsize = 22

make_maps = true
if make_maps
    for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
        if i > 1 || j > 1
            forcing = forcing_frames[i, j]

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

            scatter!(axs[1, 1], f_eq .+ f2015, V_eq, color = :red, markersize = ms2)
            text!(axs[1, 1], f_eq .+ f2015, V_eq, text = state_labels[i, j],
                color = :red, fontsize = 30, font = :bold, offset = text_offsets[i, j])

            hidedecorations!(axs[i, j])
            heatmap!(axs[i, j], xc, yc, z_bed; cmaps["z_bed2"]...)
            heatmap!(axs[i, j], xc, yc, z_srf .* f_ice; cmaps["z_srf"]...)
            contour!(axs[i, j], xc, yc, f_grnd .+ f_ice, levels = [1.9],
                color = :red, linewidth = 2)
            if xlims_frames[i, j] !== nothing
                contour!(axs[i, j], xc, yc, (xlims_frames[i, j][1] .< X .< xlims_frames[i, j][2]) .&
                    (ylims_frames[i, j][1] .< Y .< ylims_frames[i, j][2]), levels = [0.5],
                    color = :darkred, linewidth = 3)
            end
            statlab = state_labels[i, j]
            text!(axs[i, j], -2500, -2450, color = :white, font = :bold,
                text="("*statlab*")", fontsize = 30)
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
end

for k in 1:aqef.n_xps
    lines!(axs[1, 1], aqef.f[k][1:s:end] ./ polar_amplification .+ f2015, aqef.V_sle[k][1:s:end],
        linewidth = lws[k], label = xp_labels[k], color = lcolor(cycling_colors[k]))
end

scatter!(axs[1, 1], eql.f ./ polar_amplification .+ f2015, eql.V_sle;
    color = :black, label = "EQL", markersize = ms1)

axislegend(axs[1, 1], position = :lb, nbanks = 1)
relwidth = 0.8
Colorbar(fig[1, 2], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Bed elevation $z_b$ (km)", ticks = latexifyticks(-6:2:2, 1f3); cmaps["z_bed2"]...)
Colorbar(fig[1, 3], vertical = false, width = Relative(relwidth), valign = 2,
    label = L"Ice surface elevation $z_s$ (km)",
    ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); cmaps["z_srf"]...)
elem_1 = LineElement(color = :red, linewidth = 2)
elem_2 = LineElement(color = :darkred, linewidth = 2)
Legend(fig[1, 4], [elem_1, elem_2], ["Grounding line", "Highlighted region"])

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
fig
save(plotsdir("16km/hysteresis/fig1.png"), fig)
save(plotsdir("16km/hysteresis/fig1.pdf"), fig)