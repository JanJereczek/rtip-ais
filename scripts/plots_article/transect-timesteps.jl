include("../intro.jl")

pfx = datadir("output/ais/hyster/16km/aqef/")
case = "pmpt-hydro-25K/1"
file_1D = joinpath(pfx, case, "yelmo1D.nc")
file_2D = joinpath(pfx, case, "yelmo2Dsm.nc")
time_1D = ncread(file_1D, "time")
time_2D = ncread(file_2D, "time")
f = ncread(file_1D, "hyst_f_now")
i_bif_1D = findfirst(f .> 1.15)
t_bif = time_1D[i_bif]
i_bif_2D = findfirst(time_2D .> t_bif)
nx, ny = size(ncread(file_2D, "x2D"))
ncslice(file, var, i, nx, ny) = reshape(
    ncread(file, var, start = [1, 1, i], count = [-1, -1, 1]), nx, ny)
slice_fgrnd(file, i, nx, ny) = (910 / 1028) .* ncslice(file, "H_ice", i, nx, ny) .+
    ncslice(file, "z_bed", i, nx, ny) .> 0

H_ice = ncslice(file_2D, "H_ice", i_bif_2D, nx, ny)
z_bed = ncslice(file_2D, "z_bed", i_bif_2D, nx, ny)
nx, ny = size(z_bed)
x = 1:nx
y = 1:ny
X = repeat(x, 1, ny)
Y = repeat(y', nx, 1)
nxx = nyy = 150
ii = 50:50+nxx
jj = 100:100+nyy

set_theme!(theme_latexfonts())
nrows, ncols = 1, 2
fig = Figure(size=(1600, 600), fontsize = 20)
ax1 = Axis(fig[1, 1], aspect = DataAspect())
ax2 = Axis(fig[1, 2:3])
hidedecorations!(ax1)
ax2.xlabel = L"Distance along transect ($10^3$ km) $\,$"
ax2.ylabel = "Elevation (m)"
ax2.yaxisposition = :right
heatmap!(ax1, view(z_bed, ii, jj); cmaps["z_bed"]...)
heatmap!(ax1, view(z_bed .+ H_ice, ii, jj); cmaps["z_srf"]...)

t1 = 8
t2 = 18
dt = 2
t_steps = (t1:dt:t2) .* 1f3
i_steps = [argmin(abs.(time_2D .- t)) for t in t_steps]
nsnaps = length(i_steps)
grline_map = cgrad(:jet, range(0, stop = 1, length = nsnaps+1), categorical = true)
for k in eachindex(i_steps)
    contour!(ax1, 
        view(slice_fgrnd(file_2D, tsteps[k], nx, ny), ii, jj),
        levels = [0.5], color = grline_map[k], linewidth = 3)
    # scatter!(axs[2, 1], time1D[2][k * dt_scale] ./ 1e3, Vsle[2][k * dt_scale],
    #     color = grline_map[Int((k-90)/5)+1], markersize = 30)
end
Colorbar(fig[2, 1], vertical = false, width = Relative(0.45), label = "Surface elevation (km)",
    flipaxis = false, halign = :left, ticks = latexifyticks(0:4, 1f3); cmaps["z_srf"]...)
Colorbar(fig[2, 1], vertical = false, width = Relative(0.45), label = "Time (kyr)",
    flipaxis = false, halign = :right, colormap = grline_map,
    colorrange = (t_steps[1] ./ 1f3, (t_steps[end]) ./ 1f3 + dt),
    ticks = latexifyticks( Int.(t_steps ./ 1f3) ))

regions = ["ross", "abbot", "ronne"]
region = regions[2]
img = FileIO.load(datadir("processed/16km-transect-$region.png"))
rgba2gray(x) = x.alpha
img = [rgba2gray(img[i, j]) for i in 1:size(img, 1), j in 1:size(img, 2)]
mask = rotr90(img) .> 0.5
heatmap!(ax1, view(mask, ii, jj), colormap = cgrad([:black, :black]),
    colorrange = (1, 1.1), lowclip = :transparent, highclip = :transparent)
fig

dxx = 16
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
r32 = sqrt(dxx^2 + dxx^2)

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

crop_transect = 10
ordered_idx = ordered_idx[crop_transect:end-crop_transect]
r = r[crop_transect:end-crop_transect]

extract_idx(idx) = Tuple(idx.I)
ii = [extract_idx(ordered_idx[i])[1] for i in 1:length(ordered_idx)]
jj = [extract_idx(ordered_idx[i])[2] for i in 1:length(ordered_idx)]

rr = r .* 1e3

rlolim, rhilim = extrema(rr)
zlolim = -1200
zhilim = 500
xlims!(ax2, (rlolim, rhilim))
ylims!(ax2, (zlolim, zhilim))

f_grnd = slice_fgrnd(file_2D, i_bif_2D, nx, ny)
r_grline = findfirst(f_grnd[ordered_idx] .>= 1)
z_sl = zeros(nx, ny)
plot_transect(ax2, z_bed, nothing, f_grnd, H_ice, z_sl, ii, jj, rr)

# rowsize!(fig.layout, 1, 400)
# rowsize!(fig.layout, 2, 400)
# colsize!(fig.layout, 1, 400)
# colsize!(fig.layout, 2, 800)
rowgap!(fig.layout, -50)
fig
save(plotsdir("16km/bedrock/transect-viscs.png"), fig)