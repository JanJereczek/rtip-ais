include("../intro.jl")

exp1 = datadir("output/ais/spinup/base/0/yelmo2D.nc")
exp2 = datadir("output/ais/hyster/aqef/fast/0/yelmo2D.nc")
exp3 = datadir("output/ais/hyster/aqef/fast-parspinup/0/yelmo2D.nc")
# plot_diffheatmaps(exp1, exp3, "fast-parspinup", symmetric = false)

exp4 = datadir("output/ais/spinup/ensemble/long/1/yelmo2D.nc")
exp5 = datadir("output/ais/hyster/aqef/20K/1/yelmo2D.nc")
plot_diffheatmaps(exp4, exp5, "symmetric/longspinup", symmetric = true)

exp6 = datadir("output/ais/spinup/ensemble/long/1/yelmo_restart.nc")
exp7 = datadir("output/ais/hyster/aqef/20K/1/yelmo_restart_init.nc")
plot_diffheatmaps(exp6, exp7, "restarts/", symmetric = true)