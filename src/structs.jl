struct AQEFResults{T<:AbstractFloat}
    xps::Vector{String}
    n_xps::Int
    t_1D::Vector{Vector{T}}
    f::Vector{Vector{T}}
    V_sle::Vector{Vector{T}}
    A_ice::Vector{Vector{T}}
    t_2D::Vector{Vector{T}}
end

function AQEFResults(T, xps::Vector{String})
    n_xps = length(xps)
    
    t_1D = Vector{Vector{T}}(undef, n_xps)
    f = Vector{Vector{T}}(undef, n_xps)
    V_sle = Vector{Vector{T}}(undef, n_xps)
    A_ice = Vector{Vector{T}}(undef, n_xps)
    t_2D = Vector{Vector{T}}(undef, n_xps)


    for (i, xp) in enumerate(xps)
        xpdir = joinpath([xp, "0"])

        @show xpdir
        t_1D[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "time")
        f[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "hyst_f_now")
        V_sle[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "V_sle")
        A_ice[i] = ncread(joinpath(xpdir, "yelmo1D.nc"), "A_ice")
        t_2D[i] = ncread(joinpath(xpdir, "yelmo2D.nc"), "time")
    end
    AQEFResults(xps, n_xps, t_1D, f, V_sle, A_ice, t_2D)
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
    cmb = zeros(nt)
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
            cmb[k] = mb_net[k] .- bmb[k] .- smb[k]
        else
            mb_resid[k] = mean(ncslice(file, "mb_resid", k + k_offset)[mask_applied])
            cmb[k] = mean(ncslice(file, "cmb", k + k_offset)[mask_applied])
        end
    end
    mb_net = bmb .+ smb .+ cmb
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


mutable struct HeatmapRtip{T}
    visc_cases::Vector{String}
    n_visc_cases::Int
    f::Vector{Vector{T}}
    rtime::Vector{Vector{T}}
    V::Vector{Vector{T}}
    t_end::Vector{Vector{T}}
    paths::Vector{Vector{String}}
end

function HeatmapRtip{T}(visc_cases::Vector{String}) where T
    n_visc_cases = length(visc_cases)
    f = [ T[] for _ in 1:n_visc_cases ]
    rtime = [ T[] for _ in 1:n_visc_cases ]
    V = [ T[] for _ in 1:n_visc_cases ]
    t_end = [ T[] for _ in 1:n_visc_cases ]
    paths = [ String[] for _ in 1:n_visc_cases ]
    return HeatmapRtip(visc_cases, n_visc_cases, f, rtime, V, t_end, paths)
end

function aggregate_xp!(hr::HeatmapRtip, dir::String)

    (;visc_cases, f, rtime, V, t_end, paths) = hr
    params = readdlm(joinpath(dir, "info.txt"))
    i_dirs = findfirst(params[1, :] .== "runid")
    i_rtime = findfirst(params[1, :] .== "hyster.dt_ramp")
    i_f_max = findfirst(params[1, :] .== "hyster.f_max")
    i_visc_case = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
    i_scale = findfirst(params[1, :] .== "isos.viscosity_scaling")
    n_files = size(params, 1) - 1

    for j in 1:n_files
        file = joinpath(dir, "$(params[j + 1, i_dirs])", "yelmo1D.nc")
        visc_case = params[j + 1, i_visc_case]

        if visc_case == "stddev"
            if params[j + 1, i_scale] == -2.0
                visc_case = "m2stddev"
            elseif params[j + 1, i_scale] == -1.0
                visc_case = "m1stddev"
            elseif params[j + 1, i_scale] == 0.0
                visc_case = "nominal"
            elseif params[j + 1, i_scale] == 1.0
                visc_case = "p1stddev"
            elseif params[j + 1, i_scale] == 2.0
                visc_case = "p2stddev"
            end
        end
        i = findfirst(visc_cases .== visc_case)

        t_e = ncread(file, "time")[end]
        if (t_e > 10 && i !== nothing)
            push!(t_end[i], t_e)
            push!(f[i], params[j + 1, i_f_max])
            push!(rtime[i], params[j + 1, i_rtime])
            push!(V[i], ncread(file, "V_sle")[end])
            push!(paths[i], joinpath(dir, "$(params[j + 1, i_dirs])"))
        end
    end
end

function get_rtime_idx(rtime)
    rtime_vals = reverse(unique(rtime))
    rtime_idx = similar(rtime, Int)
    for i in eachindex(rtime_vals)
        inds = findall(rtime .== rtime_vals[i])
        for ind in inds
            rtime_idx[ind] = i
        end
    end
    return rtime_idx
end


mutable struct HeatmapRtipRate{T}
    visc_cases::Vector{String}
    n_visc_cases::Int
    f::Vector{Vector{T}}
    dfdt::Vector{Vector{T}}
    V::Vector{Vector{T}}
    t_end::Vector{Vector{T}}
    paths::Vector{Vector{String}}
end

function HeatmapRtipRate{T}(visc_cases::Vector{String}) where T
    n_visc_cases = length(visc_cases)
    f = [ T[] for _ in 1:n_visc_cases ]
    dfdt = [ T[] for _ in 1:n_visc_cases ]
    V = [ T[] for _ in 1:n_visc_cases ]
    t_end = [ T[] for _ in 1:n_visc_cases ]
    paths = [ String[] for _ in 1:n_visc_cases ]
    return HeatmapRtipRate(visc_cases, n_visc_cases, f, dfdt, V, t_end, paths)
end

function aggregate_xp!(hr::HeatmapRtipRate, dir::String)

    (;visc_cases, f, dfdt, V, t_end, paths) = hr
    params = readdlm(joinpath(dir, "info.txt"))
    i_dirs = findfirst(params[1, :] .== "runid")
    i_df_dt_max = findfirst(params[1, :] .== "hyster.df_dt_max")
    i_f_max = findfirst(params[1, :] .== "hyster.f_max")
    i_visc_case = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
    i_scale = findfirst(params[1, :] .== "isos.viscosity_scaling")
    n_files = size(params, 1) - 1

    for j in 1:n_files
        file = joinpath(dir, "$(params[j + 1, i_dirs])", "yelmo1D.nc")
        # @show file
        visc_case = params[j + 1, i_visc_case]

        if visc_case == "stddev"
            if params[j + 1, i_scale] == -2.0
                visc_case = "m2stddev"
            elseif params[j + 1, i_scale] == -1.0
                visc_case = "m1stddev"
            elseif params[j + 1, i_scale] == 0.0
                visc_case = "nominal"
            elseif params[j + 1, i_scale] == 1.0
                visc_case = "p1stddev"
            elseif params[j + 1, i_scale] == 2.0
                visc_case = "p2stddev"
            end
        end
        i = findfirst(visc_cases .== visc_case)

        t_e = ncread(file, "time")[end]
        if (t_e > 10 && i !== nothing)
            push!(t_end[i], t_e)
            push!(f[i], params[j + 1, i_f_max])
            push!(dfdt[i], params[j + 1, i_df_dt_max])
            push!(V[i], ncread(file, "V_sle")[end])
            push!(paths[i], joinpath(dir, "$(params[j + 1, i_dirs])"))
        end
    end
end

function get_dfdt_idx(dfdt)
    dfdt_vals = reverse(unique(dfdt))
    dfdt_idx = similar(dfdt, Int)
    for i in eachindex(dfdt_vals)
        inds = findall(dfdt .== dfdt_vals[i])
        for ind in inds
            dfdt_idx[ind] = i
        end
    end
    return dfdt_idx
end

mutable struct StepRtip{T}
    visc_cases::Vector{String}
    n_visc_cases::Int
    f::Vector{Vector{T}}
    V::Vector{Vector{T}}
    t_end::Vector{Vector{T}}
    files::Vector{Vector{String}}
end

function StepRtip{T}(visc_cases::Vector{String}) where T
    n_visc_cases = length(visc_cases)
    f = [ T[] for _ in 1:n_visc_cases ]
    V = [ T[] for _ in 1:n_visc_cases ]
    t_end = [ T[] for _ in 1:n_visc_cases ]
    files = [ String[] for _ in 1:n_visc_cases ]
    return StepRtip(visc_cases, n_visc_cases, f, V, t_end, files)
end


function aggregate_xp!(sr::StepRtip, dir)
    (;visc_cases, f, V, t_end) = sr
    params = readdlm(joinpath(dir, "info.txt"))
    i_dirs = findfirst(params[1, :] .== "rundir")
    i_f_max = findfirst(params[1, :] .== "hyster.f_max")
    i_visc_case = findfirst(params[1, :] .== "isos.viscosity_scaling_method")
    i_df_dt_max = findfirst(params[1, :] .== "hyster.df_dt_max")
    i_scale = findfirst(params[1, :] .== "isos.viscosity_scaling")
    n_files = size(params, 1) - 1

    for j in 1:n_files
        file = joinpath(dir, "$(params[j + 1, i_dirs])", "yelmo1D.nc")
        visc_case = params[j + 1, i_visc_case]

        if visc_case == "stddev"
            if params[j + 1, i_scale] == -2.0
                visc_case = "m2stddev"
            elseif params[j + 1, i_scale] == -1.0
                visc_case = "m1stddev"
            elseif params[j + 1, i_scale] == 0.0
                visc_case = "nominal"
            elseif params[j + 1, i_scale] == 1.0
                visc_case = "p1stddev"
            elseif params[j + 1, i_scale] == 2.0
                visc_case = "p2stddev"
            end
        end
        i = findfirst(visc_cases .== visc_case)

        t_e = ncread(file, "time")[end]
        if (t_e > 1_000) && (i !== nothing) && (params[j + 1, i_df_dt_max] > 0.5)
            push!(t_end[i], t_e)
            push!(f[i], params[j + 1, i_f_max] / polar_amplification + f2020)
            V[i] = vcat(V[i], ncread(file, "V_sle")[end])
            push!(sr.files[i], file)
        end
    end
end

function Base.sort!(sr::StepRtip{T}) where T
    for i in 1:sr.n_visc_cases
        idx = sortperm(sr.f[i])
        sr.f[i] = sr.f[i][idx]
        sr.V[i] = sr.V[i][idx]
        sr.t_end[i] = sr.t_end[i][idx]
        sr.files[i] = sr.files[i][idx]
    end
    return sr
end
