z_srf_map = cgrad([:gray10, :gray95])
z_bed_map = cgrad(:oleron, [2/3])
uxy_s_map = cgrad(:Purples, 0:0.25:1, categorical = true)
u_bed_map = cgrad([:orange, :white, :midnightblue], [1/9])
z_sl_map = cgrad(:BrBg)
smb_map = cgrad(:PuOr, [6/8])
T_shlf_map = cgrad(:jet)

cmaps = Dict(
    "z_srf" => (colormap = z_srf_map, colorrange = (1e1, 4e3),
        lowclip = :transparent, highclip = :white),
    "z_bed" => (colormap = z_bed_map, colorrange = (-4e3, 2e3),
        lowclip = oleronmap[1], highclip = oleronmap[end]),
    "uxy_s" => (colormap = uxy_s_map, colorrange = (0, 4),
        lowclip = :transparent, highclip = uxy_s_map[end]),
    "u_bed" => (colormap = u_bed_map, colorrange = (-100, 800),
        lowclip = u_bed_map[1], highclip = u_bed_map[end]),
    "z_sl" => (colormap = z_sl_map, colorrange = (-10, 10),
        lowclip = z_sl_map[1], highclip = z_sl_map[end]),
    "smb" => (colormap = smb_map, colorrange = (-60, 20)),
    "T_shlf" => (colormap = T_shlf_map, colorrange = (260, 280)),
    "f_grnd" => (levels = [0.5], color = :red, linewidth = 3),
)