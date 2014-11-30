module DCEMRI

using ArgParse
using Calculus    #.jacobian
using MAT

export ser, r1eff, tissueconc, fitr1, fitdce, runmodel, parsefromargs,
  defaultdict, ccc, nlsfit, rundemo, makeplots

const verbose = true
const version = "v0.3"

include("util.jl")
include("fitting.jl")
include("models.jl")
include("science.jl")
include("plotting.jl")
include("rundemo.jl")

end
