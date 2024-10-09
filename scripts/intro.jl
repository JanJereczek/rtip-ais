using DrWatson
@quickactivate "rtip-ais"

using CairoMakie
using DelimitedFiles
using Interpolations
using NCDatasets
using Statistics

# Here you may include files from the source directory
include(srcdir("dummy_src_file.jl"))

println(
"""
Currently active project is: $(projectname())

Path of active project: $(projectdir())
"""
)
