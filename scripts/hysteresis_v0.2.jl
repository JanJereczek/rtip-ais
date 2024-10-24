include("intro.jl")

nrows, ncols = 2, 4
fig = Figure(size=(1600, 980), fontsize = 24)
axs = [Axis(fig[i+1, j]) for i in 1:nrows, j in 1:ncols]

dir = datadir("output/ais/hyster/aqef/")
xps = ["20K/3", "20K/2", "regrowth/20K-3_equil-24/0"]
#, "20K-dt5/0/", "20K-dt5-reducedout/0/"] #, "20K/0/", "20K/1/", "20K/2/", "20K/3/"]

stride = 2_000
var_names_1D = ["time", "V_sl"] #, "V_ice", "A_ice"]
forcing_frames = reshape([0, 6, 6.2, 8.2, 10.5, 12], 3, 2)'
state_labels = string.(reshape(1:6, 3, 2)')

ref_xp = xps[1]
ref_file1D = get_files(joinpath(dir, ref_xp))[1]
ref_file2D = get_files(joinpath(dir, ref_xp))[5]
ref_f = load_netcdf(ref_file1D, ["hyst_f_now"], stride)[1]
ref_time = load_netcdf(ref_file1D, ["time"], stride)[1]

for l in eachindex(xps)
    xp = xps[l]
    fldr = joinpath(dir, xp)
    file1D, file1Dwais, _, file1Deais, file2D, file2Dwais, file2Dsm, _ = get_files(fldr)
    files1D = [file1Dwais, file1Deais]
    f = load_netcdf(file1D, ["hyst_f_now"], stride)[1]

    for k in eachindex(files1D)
        vars = load_netcdf(files1D[k], var_names_1D)
        time, V_sl = vars    #, V_ice, A_ice

        i1 = 1
        i2 = minimum(lastindex.(vars))

        if xp == ref_xp
            lines!(axs[k, 1], f[i1:i2], V_sl[i1:i2], label = "$l", linewidth = 4)
        
            for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
                i3 = findfirst(f .> forcing_frames[i, j])
                scatter!(axs[k, 1], f[i3], V_sl[i3], color = :red, markersize = 25)
                text!(axs[k, 1], f[i3], V_sl[i3], text = state_labels[i, j],
                    color = :grey10, fontsize = 30, font = :bold, offset = (-20, -40))
            end

        else
            lines!(axs[k, 1], f[i1:i2], V_sl[i1:i2], label = "$l", linewidth = 2)
        end
    end

    if xp == ref_xp
        T = Float32
        V_WAIS_equil = T[]
        V_EAIS_equil = T[]
        f_equil = T[]
        for (root, dirs, files) in walkdir(datadir("output/ais/hyster/equil/20K"))
            for file in files
                if occursin("yelmo1D_WAIS.nc", file)
                    push!(V_WAIS_equil, 
                        mean(ncread(joinpath(root, file), "V_sl")[end-10:end]))
                elseif occursin("yelmo1D_EAIS.nc", file)
                    push!(V_EAIS_equil, 
                        mean(ncread(joinpath(root, file), "V_sl")[end-10:end]))
                elseif occursin("yelmo1D.nc", file)
                    push!(f_equil,
                        mean(ncread(joinpath(root, file), "hyst_f_now")[1:stride:end]))
                end
            end
        end
        idx = sortperm(f_equil)
        equil_opts = (label = "equil", color = :black, linewidth = 4)
        scatterlines!(axs[1, 1], f_equil[idx], V_WAIS_equil[idx]; equil_opts...)
        scatterlines!(axs[2, 1], f_equil[idx], V_EAIS_equil[idx]; equil_opts...)
    end
end

for i in 1:nrows
    axs[i, 1].xminorticks = IntervalsBetween(5) 
    axs[i, 1].xminorgridvisible = true
    axislegend(axs[i, 1], position = :lb, fontsize = 16)
end

# axs[2, 1].xticksvisible = false
# axs[2, 1].xticklabelsvisible = false
axs[1, 1].xaxisposition = :top

axs[1, 1].xlabel = "Regional atmospheric warming (K)"
axs[nrows, 1].xlabel = "Regional atmospheric warming (K)"

# axs[1, 1].ylabel = "AIS volume (10^6 km^3)"
axs[1, 1].ylabel = "WAIS volume (10^6 km^3)"
axs[nrows, 1].ylabel = "EAIS volume (10^6 km^3)"

[hidedecorations!(axs[i, j]) for i in 1:nrows, j in 2:ncols]





_, _, _, _, file2D, file2Dwais, file2Dsm, _ = get_files(joinpath(dir, ref_xp))
var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
oleronmap = cgrad(:oleron, [2/3])
bedmap = (colormap = oleronmap, colorrange = (-4e3, 2e3),
    lowclip = oleronmap[1], highclip = oleronmap[end])
icemap = (colormap = cgrad([:gray10, :gray95]), colorrange = (1e1, 4e3), lowclip = :transparent, highclip = :white)
velmap = (colormap = :inferno, colorrange = (1e2, 1e4), lowclip = :transparent)
t2D = ncread(file2D, "time")

for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    forcing = forcing_frames[i, j]
    i3 = findfirst(ref_f .>= forcing)
    frame_index = argmin(abs.(t2D .- ref_time[i3]))

    z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D, frame_index)
    heatmap!(axs[i, j+1], z_bed; bedmap...)
    heatmap!(axs[i, j+1], z_srf .* f_ice; icemap...)
    contour!(axs[i, j+1], f_grnd .+ f_ice, levels = [1.9], color = :darkorange, linewidth = 2)
    text!(axs[i, j+1], 10, 10, color = :white, font = :bold, text=state_labels[i, j], fontsize = 30)
    # heatmap!(axs[1, 2], uxy_srf; velmap...)
end

[axs[i, j].aspect = DataAspect() for i in 1:nrows, j in 2:ncols]
relwidth = 0.4
rowgap!(fig.layout, 5)
colgap!(fig.layout, 5)
rowgap!(fig.layout, 1, -40)
Colorbar(fig[1, 2:ncols], vertical = false, width = Relative(relwidth), label = "Bed elevation (km)", ticks = (-4e3:1e3:2e3, string.(-4:2)); bedmap...)
Colorbar(fig[nrows+2, 2:ncols], vertical = false, width = Relative(relwidth), label = "Surface elevation (km)", flipaxis = false, ticks = (vcat(1e1, 1e3:1e3:4e3), string.(0:4)); icemap...)
rowgap!(fig.layout, nrows+1, -40)

fig
save(plotsdir("hysteresis/1.png"), fig)
save(plotsdir("hysteresis/1.pdf"), fig)