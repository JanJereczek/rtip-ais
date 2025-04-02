struct AQEFResults{T<:AbstractFloat}
    xps::Vector{String}
    n_xps::Int
    t_1D::Vector{Vector{T}}
    f::Vector{Vector{T}}
    V_sle::Vector{Vector{T}}
    t_2D::Vector{Vector{T}}
end

function AQEFResults(T, xps::Vector{String})
    n_xps = length(xps)
    
    t_1D = Vector{Vector{T}}(undef, n_xps)
    f = Vector{Vector{T}}(undef, n_xps)
    V_sle = Vector{Vector{T}}(undef, n_xps)

    t_2D = Vector{Vector{T}}(undef, n_xps)


    for (i, xp) in enumerate(xps)
        xpdir = joinpath([xp, "0"])

        t_1D[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "time")
        f[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "hyst_f_now")
        V_sle[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "V_sle")

        t_2D[i] = ncread(joinpath(xpdir, "yelmo2D.nc"), "time")
    end
    AQEFResults(xps, n_xps, t_1D, f, V_sle, t_2D)
end

struct MaskedMassbalance{T}
    time::Vector{T}
    bmb::Vector{T}
    smb::Vector{T}
    cmb::Vector{T}
    mb_net::Vector{T}
    mb_resid::Vector{T}
end

function masked_massbalance(file, t1, t2, mask; T = Float32)
    t = ncread(file, "time")
    k_offset = findfirst(t .>= t1) - 1
    time = t[t1 .<= t .<= t2]
    nt = length(time)
    bmb = zeros(nt)
    smb = zeros(nt)
    mb_net = zeros(nt)
    mb_resid = zeros(nt)
    for k in 1:nt
        if isnothing(mask)
            mask_applied = ncslice(file, "H_ice", k + k_offset) .> 1
        else
            mask_applied = mask
        end
        bmb[k] = mean(ncslice(file, "bmb", k + k_offset)[mask_applied])
        smb[k] = mean(ncslice(file, "smb", k + k_offset)[mask_applied])
        mb_net[k] = mean(ncslice(file, "mb_net", k + k_offset)[mask_applied])
        if occursin("sm", file)
            nothing
        else
            mb_resid[k] = mean(ncslice(file, "mb_resid", k + k_offset)[mask_applied])
        end
    end
    cmb = mb_net .- bmb .- smb
    return MaskedMassbalance(time, T.(bmb), T.(smb), T.(cmb), T.(mb_net), T.(mb_resid))
end

struct EquilResults{T<:AbstractFloat}
    dir::String
    files1D::Vector{String}
    files2D::Vector{String}
    f::Vector{T}
    V_sle::Vector{T}
    nt::Vector{Int}
end

function EquilResults(T, dir; level = 3)
    files1D = recursive_global(dir, "yelmo1D.nc", level)
    files2D = recursive_global(dir, "yelmo2D.nc", level)
    n_xps = length(files1D)
    f = zeros(T, n_xps)
    V_sle = zeros(T, n_xps)
    nt = zeros(Int, n_xps)
    for i in eachindex(files1D)
        f[i] = ncread(files1D[i], "hyst_f_now")[end]
        V_sle[i] = ncread(files1D[i], "V_sle")[end]
        nt[i] = length(ncread(files2D[i], "time"))
    end
    idx = sortperm(f)
    return EquilResults(dir, files1D[idx], files2D[idx], f[idx], V_sle[idx], nt[idx])
end

function hyst!(ax, eql::EquilResults; kwargs...)
    lines!(ax, eql.f, eql.V_sle; kwargs...)
    return nothing
end

function hm!(ax, eql::EquilResults, f; kwargs...)
    i = findfirst(eql.f .>= f)
    z_bed, z_srf, f_grnd, f_ice = load_netcdf_2D(eql.files2D[i],
        ["z_bed", "z_srf", "f_grnd", "f_ice"], eql.nt[i])
    heatmap!(ax, z_bed; kwargs...)
    heatmap!(ax, z_srf .* f_ice; kwargs...)
    contour!(ax, f_grnd .+ f_ice, levels = [1.9], color = :red, linewidth = 2)
    return nothing
end
