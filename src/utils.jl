latexify(s) = L"%$s $\,$"
latexifyticks(x) = (x, latexify.(x))
latexifyticks(x, scale) = (x .* scale, latexify.(x))

function get_files(dir)
    file1D = joinpath(dir, "yelmo1D.nc")
    file1Dwais = joinpath(dir, "yelmo1D_WAIS.nc")
    file1Dapis = joinpath(dir, "yelmo1D_APIS.nc")
    file1Deais = joinpath(dir, "yelmo1D_EAIS.nc")
    file2D = joinpath(dir, "yelmo2D.nc")
    file2Dwais = joinpath(dir, "yelmo2Dwais.nc")
    file2Dsm = joinpath(dir, "yelmo2Dsm.nc")
    file3D = joinpath(dir, "yelmo3D.nc")
    file_bsl = joinpath(dir, "bsl.nc")

    return file1D, file1Dwais, file1Dapis, file1Deais, file2D, file2Dwais, file2Dsm,
        file3D, file_bsl
end

function load_netcdf(file, vars::Vector{String})
    return [ncread(file, var) for var in vars]
end

function load_netcdf(file, vars::Vector{String}, stride::Int)
    return [ncread(file, var)[1:stride:end] for var in vars]
end

function load_netcdf_2D(file, vars::Vector{String}, index)
    return [ncread(file, var, start = [1, 1, index], count = [-1, -1, 1])[:, :, 1] for var in vars]
end

function recursive_global(dir, pattern, level)
    filepaths = String[]
    gl = "$pattern"
    for _ in 1:level
        append!(filepaths, glob(gl, dir))
        gl = "*/" * gl
    end
    return filepaths
end

# After timing this function, it turns out perf and alloc is similar to load_netcdf
function load_ncdataset_1D(file, vars)
    NCDataset(file) do ds
        return [ds[var][:] for var in vars]
    end
end


function load_ssps(scale; wrt = :pd)
    hist = readdlm(datadir("processed/SSP/History.csv"), ',')
    ssp1 = readdlm(datadir("processed/SSP/SSP2-34.csv"), ',')
    ssp2 = readdlm(datadir("processed/SSP/SSP2.csv"), ',')
    ssp3 = readdlm(datadir("processed/SSP/SSP3.csv"), ',')
    ssp5 = readdlm(datadir("processed/SSP/SSP5.csv"), ',')
    if wrt == :pd
        f2015 = hist[end, 2]
        view(ssp1, :, 2) .-= f2015
        view(ssp2, :, 2) .-= f2015
        view(ssp3, :, 2) .-= f2015
        view(ssp5, :, 2) .-= f2015
    elseif wrt == :pi
        nothing
    else
        error("wrt must be either :pd or :pi")
    end

    view(ssp1, :, 2) .*= scale
    view(ssp2, :, 2) .*= scale
    view(ssp3, :, 2) .*= scale
    view(ssp5, :, 2) .*= scale
    return ssp1, ssp2, ssp3, ssp5
end

function shade_transitions!(axs, shade, offset, alpha, color)
    for i in eachindex(shade)
        sshade = (shade[i][1]:0.01:shade[i][2]) .+ offset
        for ax in axs
            vlines!(ax, sshade, alpha = alpha, color = color)
        end
    end
end

# struct StichedData1D{T}
#     files::Vector{String}
#     t_vecs::Vector{Vector{T}}
#     t_stiched::Vector{T}
# end

# function StichedData1D(files)
#     t_vecs = [ncread(file, "time") for file in files]
#     t_stiched = vcat(t_vecs[1], t_vecs[2] .+ t_vecs[1][end])
#     return StichedData1D(files, t_vecs, t_stiched)
# end

# struct StichedData2D{T}
#     files::Vector{String}
#     t_vecs::Vector{Vector{T}}
#     t_vecs_offset::Vector{Vector{T}}
#     t_stiched::Vector{T}
#     t_end::Vector{T}
# end

# function StichedData2D(files)
#     t_vecs = [ncread(file, "time") for file in files]
#     t_vecs_offset = [t_vecs[1], t_vecs[2] .+ t_vecs[1][end]] 
#     t_stiched = vcat(t_vecs_offset...)
#     t_end = vcat(t_vecs_offset[1][end], t_vecs_offset[2][end])
#     return StichedData2D(files, t_vecs, t_vecs_offset, t_stiched, t_end)
# end

# function NetCDF.ncread(sd::StichedData1D, var::String)
#     return vcat([ncread(sd.files[i], var) for i in 1:length(sd.files)]...)
# end

# function NetCDF.ncread(sd::StichedData2D, var::String, time)
#     file_idx = findfirst(sd.t_end .>= time)
#     time_idx = findfirst(sd.t_vecs_offset[file_idx] .>= time)
#     return dropdims(ncread(sd.files[file_idx], var, start = [1, 1, time_idx],
#         count = [-1, -1, 1]), dims = 3)
# end


ncslice(file, var, idx) = dropdims(
    ncread(file, var, start = [1, 1, idx], count = [-1, -1, 1]), dims = 3)