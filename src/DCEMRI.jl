module DCEMRI

using ArgParse
using Calculus    #.jacobian
using MAT
using Pkg, LinearAlgebra, Random, Statistics, Distributed, Printf, LsqFit

export ser, r1eff, tissueconc, fitr1, fitdce, fitdata,
  defaults, ccc, nlsfit, makeplots, demo, validate

const verbose = true
const version = v"0.2.2"

include("util.jl")
include("fitting.jl")
include("models.jl")
include("science.jl")
include("plotting.jl")
include("demo.jl")
include("validate.jl")

end
