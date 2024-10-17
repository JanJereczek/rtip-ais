function get_files(dir)
    file1D = joinpath(dir, "yelmo1D.nc")
    file1Dwais = joinpath(dir, "yelmo1D_WAIS.nc")
    file1Dapis = joinpath(dir, "yelmo1D_APIS.nc")
    file1Deais = joinpath(dir, "yelmo1D_EAIS.nc")
    file2D = joinpath(dir, "yelmo2D.nc")
    file2Dwais = joinpath(dir, "yelmo2Dwais.nc")
    file2Dsm = joinpath(dir, "yelmo2Dsm.nc")
    file3D = joinpath(dir, "yelmo3D.nc")

    return file1D, file1Dwais, file1Dapis, file1Deais, file2D, file2Dwais, file2Dsm, file3D
end


function load_netcdf(file, vars::Vector{String})
    return [ncread(file, var) for var in vars]
end

function load_netcdf(file, vars::Vector{String}, stride::Int)
    return [ncread(file, var)[1:stride:end] for var in vars]
end

function load_netcdf_2D(file, vars::Vector{String}, index)
    return [ncread(file, var)[:, :, index] for var in vars]
end










# After timing this function, it turns out perf and alloc is similar to load_netcdf
function load_ncdataset_1D(file, vars)
    NCDataset(file) do ds
        return [ds[var][:] for var in vars]
    end
end