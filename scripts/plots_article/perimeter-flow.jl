include("../intro.jl")

polar_amplification = 1.8
f_to = 0.25

T = Float32
visc_type = "normvisc"    # "equil" or "aqef"
regrowth_dir = datadir("output/ais/hyster/16km/regrowth")
xp = "$regrowth_dir/aqef/pmpt-$visc_type-fastnormforcing/0"

fn1D = "$xp/yelmo1D.nc"
fn2D = "$xp/yelmo2D.nc"
f1D = ncread(fn1D, "hyst_f_now")
t1D = ncread(fn1D, "time")
t2D = ncread(fn2D, "time")
x, y = ncread(fn2D, "xc"), ncread(fn2D, "yc")
X, Y = ncread(fn2D, "x2D"), ncread(fn2D, "y2D")
data, f_ice, z_bed, z_srf = similar(X), similar(X), similar(X), similar(X)
nanmask = fill(false, size(X))

k_bif_2D = 53
# k_bif_1D = argmin(abs.(t1D .- t2D[k_bif_2D]))
t_end_bif = t1D[findfirst(f1D ./ polar_amplification .< 7.5)]
k_end_bif = argmin(abs.(t2D .- t_end_bif))
x1, x2 = 500, 1500
y1, y2 = 200, 1200

k_snaps = [k_bif_2D, k_bif_2D + 2, k_bif_2D + 4]
vars = ["taud_acy", "uy_s", "taud_acx", "ux_s", "smb", "visc_eff_int"]
varlabels = [
    L"Driving stress in $y$",
    L"Surface velocity in $y$",
    L"Driving stress in $x$",
    L"Surface velocity in $x$",
    L"Surface mass balance $\,$",
    L"Effective ice viscosity $\,$"
]
def_nmap = 20
taumap = cgrad(:PuOr, range(0, stop = 1, length = def_nmap), categorical = true)
umap = cgrad(:RdBu, range(0, stop = 1, length = def_nmap), rev = true, categorical = true)
logumap = cgrad([umap[1], umap[5], :white, :white, umap[16], umap[end]],
    range(0, stop = 1, length = def_nmap), categorical = true)
smbmap = cgrad(:BrBg, range(0, stop = 1, length = def_nmap), categorical = true)
numap = cgrad(:plasma, range(0, stop = 1, length = def_nmap), categorical = true)
umax = 500

tauopts = (colormap = taumap, colorrange = (-20, 20), lowclip = taumap[1],
    highclip = taumap[end])
uopts = (colormap = logumap, colorrange = (-log10.(umax), log10(umax)), lowclip = umap[1],
    highclip = umap[end])
smbopts = (colormap = smbmap, colorrange = (-0.5, 0.5), lowclip = smbmap[1],
    highclip = smbmap[end])
nuopts = (colormap = numap, colorrange = (10.5, 11.5), lowclip = numap[1],
    highclip = numap[end])
copts = [tauopts, uopts, tauopts, uopts, smbopts, nuopts]
scale_opts = [1, 2, 1, 2, 3, 4]

set_theme!(theme_latexfonts())
n_snaps, n_vars = length(k_snaps), length(vars)
nrows, ncols = n_snaps, n_vars
fig = Figure(size = (1600, 900), fontsize = 20)
axs = [Axis(fig[i, j], aspect = AxisAspect(1)) for i in 1:nrows, j in 1:ncols]

for i in 1:nrows
    k = k_snaps[i]
    f_ice .= ncslice(fn2D, "f_ice", k)
    z_srf .= ncslice(fn2D, "z_srf", k)
    nanmask .= (f_ice .< 0.001)
    z_srf[nanmask] .= NaN

    for j in 1:ncols

        if i == 1
            axs[i, j].title = varlabels[j]
        end

        var = vars[j]
        data .= ncslice(fn2D, var, k)
        data[nanmask] .= NaN
        
        if scale_opts[j] == 1
            heatmap!(axs[i, j], x, y, data ./ 1f3; copts[j]...)
        elseif scale_opts[j] == 2
            heatmap!(axs[i, j], x, y, log10.(abs.(data)) .* sign.(data); copts[j]...)
        elseif scale_opts[j] == 3
            heatmap!(axs[i, j], x, y, data ./ 1f0; copts[j]...)
        elseif scale_opts[j] == 4
            heatmap!(axs[i, j], x, y, log10.(data); copts[j]...)
        end

        contour!(axs[i, j], x, y, z_srf, levels = 0:200:2000, linewidth = 1, color = :gray50)
        contour!(axs[i, j], x, y, f_ice, levels = [0.1], linewidth = 2, color = :black)

        hidexdecorations!(axs[i, j])
        axs[i, j].yticksvisible = false
        axs[i, j].yticklabelsvisible = false
        axs[i, j].ygridvisible = false
        xlims!(axs[i, j], x1, x2)
        ylims!(axs[i, j], y1, y2)
    end

    axs[i, 1].ylabel = "t = $(Int(t2D[k] ./ 1f3)) kyr"
end

rw = 0.8
uticks = [-log10(umax), -log10(umax/10), -log10(umax/100), log10(umax/100), log10(umax/10), log10(umax)]
utickvals = Int.([umax, umax/10, umax/100, -umax/100, -umax/10, -umax])
uticklabels = string.(utickvals)
Colorbar(fig[nrows + 1, 1], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"$\tau_{\text{d}, y}$ (kPa)", halign = :center; tauopts...)
Colorbar(fig[nrows + 1, 2], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"$u_y$ ($\mathrm{m \, yr^{-1}}$)", halign = :center, ticks = (uticks, uticklabels),
    ; uopts...)
Colorbar(fig[nrows + 1, 3], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"$\tau_{\text{d}, x}$ (kPa)", halign = :center; tauopts...)
Colorbar(fig[nrows + 1, 4], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"$u_x$ ($\mathrm{m \, yr^{-1}}$)", halign = :center, ticks = (uticks, uticklabels),
    ; uopts...)
Colorbar(fig[nrows + 1, 5], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"$\text{SMB}$ ($\mathrm{m \, yr^{-1}}$)", halign = :center; smbopts...)
Colorbar(fig[nrows + 1, ncols], vertical = false, width = Relative(rw), flipaxis = false,
    label = L"$\mathrm{log}_{10} \, \nu_{\mathrm{eff}}$ $\mathrm{(Pa \, yr \, m)}$", halign = :center; nuopts...)

rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
fig

save(plotsdir("16km/hysteresis/perimeterflow.png"), fig)