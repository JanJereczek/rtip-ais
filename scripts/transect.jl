include("intro.jl")

function crossection(x, ii, jj)
    crossx = zeros(eltype(x), length(ii), size(x, 3))
    for i in eachindex(ii)
        crossx[i, :] = x[ii[i], jj[i], :]
    end
    return vec(crossx)
end

function plot_transect(ax, z_bed, f_grnd, H_ice, z_sl, ii, jj)
    # Preprocess fields
    mask_continent = z_bed .> 0
    mask_grounded = f_grnd .>= 1
    mask_ice = H_ice .> 0
    mask_ocean = .!(mask_continent .| mask_grounded)
    mask_shelves = mask_ocean .& mask_ice
    mask_open_ocean = mask_ocean .& .!(mask_ice)

    zsbottom = fill(0f0, size(H_ice))
    zstop = fill(0f0, size(H_ice))
    zsbottom[mask_ocean] .= (z_sl[mask_ocean] .- H_ice[mask_ocean]) .* (931 / 1028)
    zstop[mask_ocean] .= zsbottom[mask_ocean] .+ H_ice[mask_ocean]

    z_bed_cross = crossection(z_bed, ii, jj)
    H_ice_cross = crossection(H_ice, ii, jj)
    mask_grounded_cross = crossection(mask_grounded, ii, jj)
    z_sl_cross = crossection(z_sl, ii, jj)
    zsbottom_cross = crossection(zsbottom, ii, jj)
    zstop_cross = crossection(zstop, ii, jj)
    mask_open_ocean_cross = crossection(mask_open_ocean, ii, jj)

    zz = z_bed_cross + H_ice_cross .* mask_grounded_cross
    bb = z_bed_cross
    gg = z_sl_cross
    ssmask = mask_open_ocean_cross
    zzloclip = min.(zz, bb)
    zzhiclip = max.(zz, bb)
    bbloclip = min.(1.1*zlolim, bb)
    bbhiclip = max.(1.1*zlolim, bb)
    sshiclip = max.(ssmask, bb)
    zshelfbottom = zsbottom_cross
    zshelftop = zstop_cross

    lines!(ax, rr, ssmask, color = :dodgerblue4)
    band!(ax, rr, bb, sshiclip, color = :royalblue, label = L"Open ocean $\,$")
    band!(ax, rr, zzloclip, zzhiclip, color = :skyblue1, label = L"Grounded ice $\,$")
    band!(ax, rr, zshelfbottom, zshelftop, color = :lightskyblue1, label = L"Floating ice $\,$")
    lines!(ax, rr, zz, color = :deepskyblue1, linewidth = 3)

    lines!(ax, rr, z_bed_cross, color = :saddlebrown, linewidth = 1, linestyle = :dash)
    lines!(ax, rr, bb, color = :saddlebrown, linewidth = 3)
    band!(ax, rr, bbloclip, bbhiclip, color = :tan1, label = L"Bedrock $\,$")
    axislegend(ax, position = :lt, nbanks = 2, framevisible = false)
end

stride1D = 20
nrows, ncols = 2, 2
fig = Figure(size=(1400, 950), fontsize = 20)
axs = [Axis(fig[i, j]) for i in 1:nrows, j in 1:ncols]

hidedecorations!(axs[1, 1])

axs[2, 1].xlabel = L"Time (kyr) $\,$"
axs[2, 1].ylabel = L"Barystatic sea level contribution (m) $\,$"

axs[1, 2].xaxisposition = :top
axs[1, 2].yaxisposition = :right
axs[1, 2].xlabel = L"Distance along transect ($10^3$ km) $\,$"
axs[1, 2].ylabel = L"Elevation (m) $\,$"
axs[2, 2].yaxisposition = :right
axs[2, 2].xlabel = L"Distance along transect ($10^3$ km) $\,$"
axs[2, 2].ylabel = L"Elevation (m) $\,$"

fig

dir = datadir("output/ais/ramps/wais-sparse/")
files1D = [
    joinpath(dir, "0/yelmo1D_WAIS.nc"),
    joinpath(dir, "1/yelmo1D_WAIS.nc"),
]
files2D = [
    joinpath(dir, "0/yelmo2Dwais.nc"),
    joinpath(dir, "1/yelmo2Dwais.nc"),
]


time1D = [ncread(file1D, "time")[1:stride1D:end] for file1D in files1D]
time2D = [ncread(file2D, "time") for file2D in files2D]
Vsle = [ncread(file1D, "V_sle")[1:stride1D:end] for file1D in files1D]
Vsle = [Vsl[1] .- Vsl for Vsl in Vsle]
lines!(axs[2, 1], time1D[1] ./ 1e3, Vsle[1], label = "0")
lines!(axs[2, 1], time1D[2] ./ 1e3, Vsle[2], label = "1")
axislegend(axs[2, 1], position = :lt)

t1 = 3e3
i1_1D = findfirst(time1D[1] .>= t1)
i1_2D = findfirst(time2D[1] .>= t1)

H_ice = ncread(files2D[1], "H_ice")[:, :, i1_2D]
z_bed = ncread(files2D[1], "z_bed")[:, :, i1_2D]
z_sl = ncread(files2D[1], "z_sl")[:, :, i1_2D]
f_grnd = ncread(files2D[1], "f_grnd")[:, :, i1_2D]

H_ice_2 = ncread(files2D[2], "H_ice")[:, :, i1_2D]
z_bed_2 = ncread(files2D[2], "z_bed")[:, :, i1_2D]
z_sl_2 = ncread(files2D[2], "z_sl")[:, :, i1_2D]
f_grnd_2 = ncread(files2D[2], "f_grnd")[:, :, i1_2D]

transect_aspect = 2
axs[1, 1].aspect = DataAspect()
axs[2, 1].aspect = AxisAspect(1)
axs[1, 2].aspect = AxisAspect(transect_aspect)
axs[2, 2].aspect = AxisAspect(transect_aspect)
heatmap!(axs[1, 1], z_bed; cmaps["z_bed"]...)
heatmap!(axs[1, 1], z_bed .+ H_ice; cmaps["z_srf"]...)

nx, ny = size(z_bed)
x = 1:nx
y = 1:ny
X = repeat(x, 1, ny)
Y = repeat(y', nx, 1)

# Define transect with 2 points
xp = [15, 40]
yp = [25, 48]
scatterlines!(axs[1, 1], xp, yp, color = :red, linewidth = 3, markersize = 25)

XX = hcat(xp, ones(length(xp)))'
M = inv(XX * XX' ) * XX
m, p = M * yp
xtransect = collect(xp[1]:xp[2])
ytransect = m * xtransect .+ p

ii = zeros(Int, length(xtransect))
jj = zeros(Int, length(xtransect))
cartesian_transect = CartesianIndex[]

for i in eachindex(xtransect)
    cindex = argmin( (xtransect[i] .- X).^2 + (ytransect[i] .- Y).^2 )
    iijj = Tuple(cindex)
    push!(cartesian_transect, cindex)
    ii[i] = iijj[1]
    jj[i] = iijj[2]
end

rr = sqrt.((x[ii] .- x[ii[1]]) .^ 2 + (y[jj] .- y[jj[1]]) .^ 2) .* 32e3
rlolim, rhilim = extrema(rr)
zlolim = minimum(z_bed_cross) - 1e2
zhilim = maximum(z_bed_cross .+ H_ice_cross) + 1e2

xlims!(axs[1, 2], (rlolim, rhilim))
ylims!(axs[1, 2], (zlolim, zhilim))
xlims!(axs[2, 2], (rlolim, rhilim))
ylims!(axs[2, 2], (zlolim, zhilim))

plot_transect(axs[1, 2], z_bed, f_grnd, H_ice, z_sl, ii, jj)
plot_transect(axs[2, 2], z_bed_2, f_grnd_2, H_ice_2, z_sl_2, ii, jj)

rowsize!(fig.layout, 1, 400)
rowsize!(fig.layout, 2, 400)
colsize!(fig.layout, 1, 400)
colsize!(fig.layout, 2, 800)
save(plotsdir("transect.png"), fig)