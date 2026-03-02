case = "rsb"

dirs = ["data/output/ais/v2/ramps/rsb/3/$k" for k in [25, 27, 24]]
labels = [L"$\eta - 2 \, \sigma$", L"$\eta$", L"$\eta + 2 \, \sigma$"]

file_ref = joinpath(dirs[2], "yelmo2Dsm.nc")
x = ncread(file_ref, "xc")
y = ncread(file_ref, "yc")
i1, i2, i3, i4 = [120, 170, 195, 230]
j1, j2, j3, j4 = [220, 252, 272, 340]
x1, x2, x3, x4 = x[i1], x[i2], x[i3], x[i4]
y1, y2, y3, y4 = y[j1], y[j2], y[j3], y[j4]
j23 = (j2+j3) ÷ 2
y23 = y[j23]
dx2 = 5
lk = [5, 5, 1]
# J23 = (j2+j3) ÷ 2
# Y23 = y[J23]
# for i in -10:10
#     @show
    # j23 = J23 + i
    # y23 = y[j23]
# for dk in 0:4
    # @show dk
    # dk = 0
    set_theme!(theme_latexfonts())
    fig_tra = Figure(size = (1600, 1630), fontsize = 44)
    nrows, ncols = 3, 2
    axs_tr = [Axis(fig_tra[i+1, 1], aspect = AxisAspect(3)) for i in 1:nrows]
    axs_t = [Axis(fig_tra[i+1, 1], aspect = AxisAspect(3)) for i in 1:nrows]
    rw = 0.6
    rh = 0.6
    k = 12

    for i in eachindex(dirs)
        file_2D = joinpath(dirs[i], "yelmo2Dsm.nc")
        file_bsl = joinpath(dirs[i], "bsl.nc")
        H_ice = ncread(file_2D, "H_ice", start = [1, 1, k], count = [-1, -1, 1])[:, :, 1]
        z_bed = ncread(file_2D, "z_bed", start = [1, 1, k], count = [-1, -1, 1])[:, :, 1]
        # heatmap!(axs_hm[i], x, y, z_bed; cmaps["z_bed"]...)
        # heatmap!(axs_hm[i], x, y, H_ice + z_bed; cmaps["z_srf"]...)
        # xlims!(axs_hm[i], (x1, x4))
        # ylims!(axs_hm[i], (y1, y4))
        # hidedecorations!(axs_hm[i])
        # lines!(axs_hm[i], [x2, x3], [y23, y23], color = :red, linewidth = 8)
        # axs_hm[i].ylabel = labels[i]

        H_ice_vec = ncread(file_2D, "H_ice", start = [i2, j23, k], count = [i3-i2+1, 1, 1])[:, 1, 1]
        z_bed_vec = ncread(file_2D, "z_bed", start = [i2, j23, k], count = [i3-i2+1, 1, 1])[:, 1, 1]
        z_srf_vec = z_bed_vec .+ H_ice_vec

        time2D = ncread(file_2D, "time", start = [k], count = [1])[1]
        time_bsl = ncread(file_bsl, "time")
        k_bsl = findfirst(time_bsl .>= time2D)
        bsl_snapshot = ncread(file_bsl, "bsl", start = [k_bsl], count = [1])[1]
        bslvec = fill(bsl_snapshot, length(H_ice_vec))
        # text!(axs_hm[i], x1 + 50, y4 - 320, text = "t = $(Int(round(time2D / 1e3, digits=1))) kyr", color = :white)

        z_bed_0_vec = ncread(file_ref, "z_bed", start = [i2, j23, 1], count = [i3-i2+1, 1, 1])[:, 1, 1]
        band!(axs_tr[i], x[i2:i3], z_bed_vec, bslvec, color = :royalblue)
        band!(axs_tr[i], x[i2:i3], z_bed_vec, z_srf_vec, color = :gray90)
        # lines!(axs_tr[i], x[i2:i3], z_srf_vec, color = :gray80, linewidth = 4, label = "ice surface")

        band!(axs_tr[i], x[i2:i3], min.(z_bed_vec, -1000), z_bed_vec, color = :tan1)
        lines!(axs_tr[i], x[i2:i3], z_bed_vec, color = :saddlebrown, linewidth = 4, label = "bedrock")
        lines!(axs_tr[i], x[i2:i3], z_bed_0_vec, color = :saddlebrown, linewidth = 4, linestyle = :dash, label = "initial bed")
        hlines!(axs_tr[i], bsl_snapshot, color = :midnightblue, linewidth = 4, linestyle = :dash, label = "sea level")
        xlims!(axs_tr[i], (x2, x3))
        ylims!(axs_tr[i], (-800, 800))
        axs_tr[i].yaxisposition = :right
        axs_tr[i].yticks = -600:200:800
        axs_tr[i].ylabel = "Elevation (m)"
        axs_tr[i].xgridvisible = false
        axs_tr[i].ygridvisible = false
        if i < nrows
            axs_tr[i].xticklabelsvisible = false
        else
            axs_tr[i].xlabel = L"$x$ (km)"
        end

        hidedecorations!(axs_t[i])
        xlims!(axs_t[i], (x2, x3))
        ylims!(axs_t[i], (-800, 800))
        grays = [:gray80, :gray65, :gray50, :gray35, :gray20]
        for dk in 0:4
            H_ice_vec = ncread(file_2D, "H_ice", start = [i2, j23, k+lk[i]*dk], count = [i3-i2+1, 1, 1])[:, 1, 1]
            z_bed_vec = ncread(file_2D, "z_bed", start = [i2, j23, k+lk[i]*dk], count = [i3-i2+1, 1, 1])[:, 1, 1]
            z_srf_vec = z_bed_vec .+ H_ice_vec
            l = findfirst(H_ice_vec .> 0)
            z_srf_vec[1:l-2] .= NaN
            lines!(axs_t[i], x[i2:i3], z_srf_vec, color = grays[dk+1], linewidth = 4, label = "t = $(Int(round((time2D + lk[i]*dk*1e3) / 1e3, digits=1))) kyr")
        end
        axislegend(axs_t[i], "Ice surface", position = :lt, labelsize = 24, titlesize = 24)

        # text!(axs_tr[i], x2 + dx2, 200, text = "t = $(Int(round(time2D / 1e3, digits=1))) kyr", color = :black)
        text!(axs_tr[i], x2 + dx2, -800, text = "A", color = :red, font = :bold)
        text!(axs_tr[i], x3 - 15, -800, text = "B", color = :red, font = :bold)
    end
    # Colorbar(fig_tra[1, 1], label = "Bed elevation (km)", vertical = false, width = Relative(rw), height = Relative(rh), ticks = (-4000:1000:3000, string.(-4:3)); cmaps["z_bed"]...)
    # Colorbar(fig_tra[1, 2], label = "Ice surface elevation (km)", vertical = false, width = Relative(rw), height = Relative(rh), ticks = (-4000:1000:3000, string.(-4:3)); cmaps["z_srf"]...)
    Legend(fig_tra[1, 1], axs_tr[1], nbanks = 4, valign = :bottom, patchsize = (40f0, 20f0))

    # text!(axs_hm[1], x1 + 50, y4 - 200, text = L"(a) $-2 \, \sigma$", color = :white, font = :bold)
    # text!(axs_hm[2], x1 + 50, y4 - 200, text = L"(c) $0 \, \sigma$", color = :white, font = :bold)
    # text!(axs_hm[3], x1 + 50, y4 - 200, text = L"(e) $+2 \, \sigma$", color = :white, font = :bold)
    text!(axs_tr[1], x2 + dx2, -200, text = L"a $-2 \, \sigma$", color = :white, font = :bold)
    text!(axs_tr[2], x2 + dx2, -200, text = L"b $0 \, \sigma$", color = :white, font = :bold)
    text!(axs_tr[3], x2 + dx2, -200, text = L"c $+2 \, \sigma$", color = :white, font = :bold)

    colsize!(fig_tra.layout, 1, 1350)
    rowsize!(fig_tra.layout, 1, 60)
    rowgap!(fig_tra.layout, 0)
    fig_tra

    save(plotsdir("v2/rtip/transect-rsb.png"), fig_tra)
# end