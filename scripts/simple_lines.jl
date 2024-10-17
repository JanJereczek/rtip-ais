include("intro.jl")

linesdir = plotsdir("lines")
isdir(linesdir) || mkdir(linesdir)

xp = datadir("output/ais/hyster/aqef/fast-parspinup/0/yelmo1D.nc")
# xp = datadir("output/ais/ctrl/spinup/base/0/yelmo1D.nc")
target_dir = "$linesdir/fast-parspinup"
isdir(target_dir) || mkdir(target_dir)
stride = 2_000

plot_lines(xp, target_dir, stride = stride)