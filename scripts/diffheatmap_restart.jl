include("intro.jl")

exp1 = datadir("output/ais/spinup/base/0/yelmo2D.nc")
exp2 = datadir("output/ais/hyster/aqef/fast/0/yelmo2D.nc")
exp3 = datadir("output/ais/hyster/aqef/fast-parspinup/0/yelmo2D.nc")
plot_diffheatmaps(exp1, exp3, "fast-parspinup")