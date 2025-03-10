include("../intro.jl")

function plot_hysteresis(retreat_key)
    set_theme!(theme_latexfonts())
    heatmap_frames = "aqef"    # "equil" or "aqef"
    f_to = 0.25
    nrows, ncols = 3, 4
    fig = Figure(size=(1700, 1320), fontsize = 24)
    axs = [Axis(fig[i+1, j], aspect = AxisAspect(1)) for i in 1:nrows, j in 1:ncols]

    dir = datadir("output/ais/hyster/16km/retreat/aqef/$retreat_key")
    xp = "0"
    xpdir = joinpath(dir, xp)
    @show xpdir

    var_names_1D = ["time", "V_sle"]
    forcing_frames = reshape(collect(0:0.1:1.1), 4, 3)'
    # forcing_frames = reshape([0, 0, 0.2, 0.3, 1.7, 1.8, 2.5, 2.6, 2.9, 3.2, 3.4, 3.5], 4, 3)'
    # forcing_frames = reshape([0, 0, 0.25, 0.35, 1.9, 2.0, 2.6, 2.7, 3.1, 3.2, 3.5, 3.6], 4, 3)'
    state_labels = latexify.(reshape(0:11, 4, 3)')

    file1D, file1Dwais, file1Dapis, file1Deais, file2D, file2Dwais, file2Dsm, _, file_bsl = get_files(xpdir)
    files1D = [file1Dapis, file1Dwais, file1Deais]
    f = ncread(file1D, "hyst_f_now") .* f_to
    t_1D = ncread(file1D, "time")
    t_end = maximum(t_1D)

    nx, ny = size(ncread(file2D, "x2D"))
    crop = 20
    ii = crop+1:nx-crop
    jj = crop+1:ny-crop

    k = 1
    s = 50
    time, V_sle = load_netcdf(file1D, var_names_1D)
    i2 = lastindex(t_1D)-1 # minimum(lastindex.(vars))
    lines!(axs[k, 1], f[1:s:i2], V_sle[1:s:i2], linewidth = 4, label = "AQEF")

    file1D_regrowth = datadir("output/ais/hyster/16km/regrowth/aqef/0/yelmo1D.nc")
    time_regrowth, V_sle_regrowth = load_netcdf(file1D_regrowth, var_names_1D)
    f_regrowth = ncread(file1D_regrowth, "hyst_f_now") .* f_to
    lines!(axs[k, 1], f_regrowth, V_sle_regrowth, linewidth = 4, color = Cycled(1))

    # dir_garbe = datadir("processed/garbe2020/")
    # files_garbe = readdir(dir_garbe)
    # garbe_f_to = 0.7 / 1.8
    # for j in eachindex(files_garbe)
    #     file_garbe = joinpath(dir_garbe, files_garbe[j])
    #     data, head = readdlm(file_garbe, ',', header = true)
    #     f_garbe = view(data, :, 1)
    #     V_garbe = view(data, :, 2)
    #     if occursin("equil", file_garbe)
    #         scatterlines!(axs[k, 1], f_garbe .* garbe_f_to, V_garbe, linewidth = 4, color = Cycled(3))
    #     else
    #         scatterlines!(axs[k, 1], f_garbe .* garbe_f_to, V_garbe, linewidth = 4, color = Cycled(4))
    #     end
    # end

    for i in 1:1
        axs[i, 1].xminorticks = IntervalsBetween(5)
        axs[i, 1].yminorticks = IntervalsBetween(10)
        axs[i, 1].xminorgridvisible = true
        axs[i, 1].yminorgridvisible = true
        axs[i, 1].xticks = latexifyticks(0:5)
        # axislegend(axs[i, 1], position = :lb, fontsize = 16)
    end
    axs[1, 1].xaxisposition = :top
    axs[1, 1].yticks = latexifyticks(0:10:60)
    axs[1, 1].xlabel = L"Regional oceanic warming (K) $\,$"
    axs[1, 1].ylabel = L"AIS $V_\mathrm{af}$ (m SLE)"
    ylims!(axs[1, 1], 0, 60)

    T = Float32
    V_equil = T[]
    f_equil = T[]
    file_2D_equil = String[]
    for (root, dirs, files) in walkdir(datadir("output/ais/hyster/16km/equil/25K-100"))
        for file in files
            if occursin("yelmo1D.nc", file)
                push!(V_equil, 
                    mean(ncread(joinpath(root, file), "V_sle")[end-10:end]))
                push!(f_equil,
                    mean(ncread(joinpath(root, file), "hyst_f_now")) .* f_to)
            elseif occursin("yelmo2D.nc", file)
                push!(file_2D_equil, joinpath(root, file))
            end
        end
    end

    idx = sortperm(f_equil)
    f_equil = f_equil[idx]
    V_equil = V_equil[idx]
    file_2D_equil = file_2D_equil[idx]
    equil_opts = (label = "Equilibrium", color = :black, linewidth = 4)
    scatterlines!(axs[1, 1], f_equil, V_equil; equil_opts...)
    ms = 15

    var_names_2D = ["z_bed", "z_srf", "uxy_s", "f_grnd", "f_ice"]
    oleronmap = cgrad(:oleron, [3/4])
    bedmap = (colormap = oleronmap, colorrange = (-6e3, 2e3),
        lowclip = oleronmap[1], highclip = oleronmap[end])
    icemap = (colormap = cgrad([:gray10, :gray95], 0:0.05:1, categorical =true),
        colorrange = (1e-1, 4e3), lowclip = :transparent, highclip = :white)
    velmap = (colormap = :inferno, colorrange = (1e2, 1e4), lowclip = :transparent)
    t2D = ncread(file2D, "time")

    # for i in axes(forcing_frames, 1), j in axes(forcing_frames, 2)
    #     if i > 1 || j > 1
    #         forcing = forcing_frames[i, j]
    #         if heatmap_frames == "aqef" || (i == 1 && j == 2)
    #             i3 = findfirst(f .>= forcing)
    #             f_eq, V_eq = f[i3], V_sle[i3]
    #             frame_index = argmin(abs.(t2D .- t_1D[i3]))
    #             z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file2D, var_names_2D, frame_index)
    #         else
    #             i3 = argmin((f_equil .- forcing) .^ 2)
    #             equil_dir = joinpath(datadir("output/ais/hyster/16km/equil/25K-100"),
    #                 string(i3), "0")
    #             f_eq, V_eq = f_equil[i3], V_equil[i3]
    #             frame_index = length(ncread(file_2D_equil[i3], "time"))
    #             z_bed, z_srf, uxy_srf, f_grnd, f_ice = load_netcdf_2D(file_2D_equil[i3],
    #                 var_names_2D, frame_index)
    #         end
    #         @show i3, f_eq, V_eq

    #         scatter!(axs[k, 1], f_eq, V_eq, color = :red, markersize = ms)
    #         text!(axs[k, 1], f_eq, V_eq, text = state_labels[i, j],
    #             color = :grey10, fontsize = 30, font = :bold, offset = (-20, -40))
            
    #         if mod(j, 2) == 1
    #             vlines!(axs[1, 1], forcing:0.01:forcing_frames[i, j+1], alpha = 0.3, color = :gray)
    #         end

    #         hidedecorations!(axs[i, j])
    #         heatmap!(axs[i, j], view(z_bed, ii, jj); bedmap...)
    #         heatmap!(axs[i, j], view(z_srf .* f_ice, ii, jj); icemap...)
    #         contour!(axs[i, j], view(f_grnd .+ f_ice, ii, jj), levels = [1.9], color = :red, linewidth = 2)
    #         text!(axs[i, j], 10, 10, color = :white, font = :bold, text=state_labels[i, j], fontsize = 30)
    #     end
    #     # heatmap!(axs[1, 2], uxy_srf; velmap...)
    # end

    hist = readdlm(datadir("processed/SSP/History.csv"), ',')
    f2014 = hist[end, 2]
    ssp1 = readdlm(datadir("processed/SSP/SSP1.csv"), ',') .- f2014
    ssp2 = readdlm(datadir("processed/SSP/SSP2.csv"), ',') .- f2014
    ssp3 = readdlm(datadir("processed/SSP/SSP3.csv"), ',') .- f2014
    # ssp4 = readdlm(datadir("processed/SSP/SSP4.csv"), ',')
    ssp5 = readdlm(datadir("processed/SSP/SSP5.csv"), ',') .- f2014

    polar_amplification = 1.8
    ssp1_2100 = ssp1[end, 2] * polar_amplification * f_to
    ssp2_2100 = ssp2[end, 2] * polar_amplification * f_to
    ssp3_2100 = ssp3[end, 2] * polar_amplification * f_to
    ssp5_2100 = ssp5[end, 2] * polar_amplification * f_to

    line_opts = (linewidth = 3, linestyle = :dash)
    vlines!(axs[1, 1], [ssp1_2100], color = :darkblue, label = "SSP1-2100"; line_opts...)
    vlines!(axs[1, 1], [ssp2_2100], color = :lightblue, label = "SSP2-2100"; line_opts...)
    vlines!(axs[1, 1], [ssp3_2100], color = :orange, label = "SSP3-2100"; line_opts...)
    vlines!(axs[1, 1], [ssp5_2100], color = :darkred, label = "SSP5-2100"; line_opts...)

    axislegend(axs[1, 1], position = :rt, labelsize = 20)
    relwidth = 0.7
    Colorbar(fig[1, 2:3], vertical = false, width = Relative(relwidth), halign = :left,
        label = L"Bed elevation (km) $\,$", ticks = latexifyticks(-6:2, 1f3); bedmap...)
    Colorbar(fig[1, 3:4], vertical = false, width = Relative(relwidth),
        label = L"Ice surface elevation (km) $\,$", halign = :right,
        ticks = (vcat([1], 1000:1000:4000), latexify.(0:4)); icemap...)

    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 1, -50)
    # rowgap!(fig.layout, nrows+1, -50)

    fig

    save(plotsdir("16km/hysteresis/$retreat_key-$heatmap_frames.png"), fig)
    # save(plotsdir("16km/hysteresis/$heatmap_frames.pdf"), fig)
end

retreat_keys = [
    "pmpt-lowvisc-atmforcing",
    "pmpt-lowvisc-garbeforcing",
    "pmpt-lowvisc-normforcing",
    "pmpt-lowvisc-normforcing-withrestarts",
    "pmpt-lowvisc-ocnforcing",
    "pmpt-normvisc-atmforcing",
    "pmpt-normvisc-fastnormforcing",
    "pmpt-normvisc-fastnormforcing-withrestarts",
    "pmpt-normvisc-garbeforcing",
    "pmpt-normvisc-normforcing",
    "pmpt-normvisc-normforcing-withrestarts",
    "pmpt-normvisc-ocnforcing",
]

for retreat_key in retreat_keys
    plot_hysteresis(retreat_key)
end