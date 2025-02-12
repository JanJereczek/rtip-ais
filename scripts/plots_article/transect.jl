include("../intro.jl")

function crossection(x, ii, jj)
    crossx = zeros(eltype(x), length(ii), size(x, 3))
    for i in eachindex(ii)
        crossx[i, :] = x[ii[i], jj[i], :]
    end
    return vec(crossx)
end

function plot_transect(ax, z_bed, z_bed_ref, f_grnd, H_ice, z_sl, ii, jj, rr)
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
    if isnothing(z_bed_ref)
        println("not plotting ref bedrock")
    else
        z_bed_2_cross = crossection(z_bed_ref, ii, jj)
        lines!(ax, rr, z_bed_2_cross, color = :saddlebrown, linewidth = 3, linestyle = :dash)
    end
    axislegend(ax, position = :lt, nbanks = 1, framevisible = false)
    vlines!(ax, [rr[r_grline]], color = :black, linestyle = :dash, linewidth = 4)
    rsl = round(z_sl_cross[r_grline] - z_bed_cross[r_grline], digits = 1)
    text!(ax, rr[r_grline], z_bed_cross[r_grline], text = "Relative sea level = $rsl m",
        offset = (-300, -100), font = :bold)
end

set_theme!(theme_latexfonts())
stride1D = 20
time_split = 4.2
n_1D = 1500
nrows, ncols = 2, 2
fig = Figure(size=(1400, 950), fontsize = 20)
axs = [Axis(fig[i, j]) for i in 1:nrows, j in 1:ncols]

hidedecorations!(axs[1, 1])

axs[2, 1].xlabel = "Time (kyr)"
axs[2, 1].ylabel = L"WAIS $V_\mathrm{af}$ (m SLE)"
axs[1, 2].xaxisposition = :top
axs[1, 2].yaxisposition = :right
axs[1, 2].xlabel = L"Distance along transect ($10^3$ km) $\,$"
axs[1, 2].ylabel = "Elevation (m)"
axs[2, 2].yaxisposition = :right
axs[2, 2].xlabel = L"Distance along transect ($10^3$ km) $\,$"
axs[2, 2].ylabel = "Elevation (m)"

# pfx = datadir("output/ais/ramps/")
# cases = ["wais-sparse-nomvisc/12", "wais-sparse-hivisc/10"]
pfx = datadir("output/ais/hyster/16km/aqef/")
cases = ["pmpt-hydro-25K/1", "pmpt-hydro-25K/0"]
files1D = [joinpath(pfx, case, "yelmo1D.nc") for case in cases]
files2D = [joinpath(pfx, case, "yelmo2Dsm.nc") for case in cases]

lw = 4
time1D = [ncread(file1D, "time")[1:n_1D] for file1D in files1D]
time2D = [ncread(file2D, "time") for file2D in files2D]
Vsle = [ncread(file1D, "V_sle")[1:n_1D] for file1D in files1D]
lines!(axs[2, 1], time1D[1] ./ 1e3, Vsle[1], label = "nominal viscosity", linewidth = lw)
lines!(axs[2, 1], time1D[2] ./ 1e3, Vsle[2], label = "high viscosity", color = :red, linewidth = lw)
axislegend(axs[2, 1], position = :lb)
xlims!(axs[2, 1], (0, 7))
vlines!(axs[2, 1], [time_split], color = :black, linestyle = :dash, linewidth = lw)

t1 = time_split * 1e3
i1_1D = findfirst(time1D[1] .>= t1)
i1_2D = findfirst(time2D[1] .>= t1)
nx, ny = size(ncread(files2D[1], "x2D"))

ncslice(file, var, i, nx, ny) = reshape(
    ncread(file, var, start = [1, 1, i], count = [-1, -1, 1]), nx, ny)
slice_fgrnd(file, i, nx, ny) = (910 / 1028) .* ncslice(file, "H_ice", i, nx, ny) .+
    ncslice(file, "z_bed", i, nx, ny) .> 0
H_ice = ncslice(files2D[1], "H_ice", i1_2D, nx, ny)
z_bed = ncslice(files2D[1], "z_bed", i1_2D, nx, ny)

# fgrnd_test = slice_fgrnd(files2D[1], i1_2D, nx, ny)
# fgrnd_test2 = slice_fgrnd(files2D[2], 700, nx, ny)
# heatmap(fgrnd_test)
# heatmap(fgrnd_test2)
nxx = nyy = 150
ii = 50:50+nxx
jj = 100:100+nyy
heatmap!(axs[1, 1], view(z_bed, ii, jj); cmaps["z_bed"]...)
heatmap!(axs[1, 1], view(z_bed .+ H_ice, ii, jj); cmaps["z_srf"]...)
fig

Colorbar(fig[0, 1], vertical = false, width = Relative(0.8), label = "Surface elevation";
    cmaps["z_srf"]...)
rowgap!(fig.layout, 1, -50)
nsnaps = 10
grline_map = cgrad(:jet, range(0, stop = 1, length = nsnaps), categorical = true)
dt_scale = 10
tsteps = range(45, step = 5, length = nsnaps - 1)
for k in eachindex(tsteps)
    contour!(axs[1, 1], 
        view(slice_fgrnd(files2D[1], tsteps[k], nx, ny), ii, jj),
        levels = [0.5], color = grline_map[k], linewidth = 2)
    # scatter!(axs[2, 1], time1D[2][k * dt_scale] ./ 1e3, Vsle[2][k * dt_scale],
    #     color = grline_map[Int((k-90)/5)+1], markersize = 30)
end

transect_aspect = 2
axs[1, 1].aspect = DataAspect()
axs[2, 1].aspect = AxisAspect(1)
axs[1, 2].aspect = AxisAspect(transect_aspect)
axs[2, 2].aspect = AxisAspect(transect_aspect)
fig

regions = ["ross", "abbot", "ronne"]
region = regions[3]
img = FileIO.load(datadir("16km-transect-$region.png"))
rgba2gray(x) = x.alpha
img = [rgba2gray(img[i, j]) for i in 1:size(img, 1), j in 1:size(img, 2)]
mask = rotr90(img) .> 0.5
heatmap!(axs[1, 1], view(mask, ii, jj), colormap = cgrad([:black, :black]),
    colorrange = (1, 1.1), lowclip = :transparent, highclip = :transparent)
fig

H_ice_2 = ncread(files2D[2], "H_ice")[:, :, i1_2D]
z_bed_2 = ncread(files2D[2], "z_bed")[:, :, i1_2D]
z_sl_2 = ncread(files2D[2], "z_sl")[:, :, i1_2D]
f_grnd_2 = ncread(files2D[2], "f_grnd")[:, :, i1_2D]
f_grnd_3D = ncread(files2D[2], "f_grnd")




nx, ny = size(z_bed)
x = 1:nx
y = 1:ny
X = repeat(x, 1, ny)
Y = repeat(y', nx, 1)



idx = CartesianIndices(mask)
mask_idx = idx[mask]
ordered_idx = similar(mask_idx)
r = zeros(length(mask_idx))

ordered_idx[1] = mask_idx[1]
ip1 = CartesianIndex(1, 0)
im1 = CartesianIndex(-1, 0)
jp1 = CartesianIndex(0, 1)
jm1 = CartesianIndex(0, -1)
ip1jp1 = CartesianIndex(1, 1)
im1jp1 = CartesianIndex(-1, 1)
r32 = sqrt(32^2 + 32^2)

for i in eachindex(ordered_idx)[1:end-1]
    continue_search = true
    for j in eachindex(mask_idx)
        if ((mask_idx[j] == mask_idx[i] + ip1jp1) ||
            (mask_idx[j] == mask_idx[i] + im1jp1)) && continue_search
            ordered_idx[i + 1] = mask_idx[j]
            # @show i + 1, r32
            r[i + 1] = r[i] + r32
            continue_search = false
        end
    end

    if continue_search
        for j in eachindex(mask_idx)
            if ((mask_idx[j] == mask_idx[i] + ip1) ||
                (mask_idx[j] == mask_idx[i] + im1) ||
                (mask_idx[j] == mask_idx[i])) && ((mask_idx[j] == mask_idx[i]) ||
                (mask_idx[j] == mask_idx[i] + jp1)) && continue_search
                ordered_idx[i + 1] = mask_idx[j]
                # @show i + 1, 32
                r[i + 1] = r[i] + 32
                continue_search = false
            end
        end
    end
end

extract_idx(idx) = Tuple(idx.I)
ii = [extract_idx(ordered_idx[i])[1] for i in 1:length(ordered_idx)]
jj = [extract_idx(ordered_idx[i])[2] for i in 1:length(ordered_idx)]

rr = r .* 1e3
r_grline = findfirst(rr .> 5.6e5)

rlolim, rhilim = extrema(rr)
zlolim = -2200
zhilim = 1600

xlims!(axs[1, 2], (rlolim, rhilim))
ylims!(axs[1, 2], (zlolim, zhilim))
xlims!(axs[2, 2], (rlolim, rhilim))
ylims!(axs[2, 2], (zlolim, zhilim))

plot_transect(axs[1, 2], z_bed, z_bed_2, f_grnd, H_ice, z_sl, ii, jj, rr)
plot_transect(axs[2, 2], z_bed_2, nothing, f_grnd_2, H_ice_2, z_sl_2, ii, jj, rr)

rowsize!(fig.layout, 1, 400)
rowsize!(fig.layout, 2, 400)
colsize!(fig.layout, 1, 400)
colsize!(fig.layout, 2, 800)
save(plotsdir("transect-viscs.png"), fig)

#=
# Define transect with 2 points
xp = [15, 40]
yp = [25, 48]
# scatterlines!(axs[1, 1], xp, yp, color = :red, linewidth = 3, markersize = 25)

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
=#