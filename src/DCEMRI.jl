module DCEMRI

using ArgParse
using Calculus    #.jacobian
using MAT

export ser, r1eff, tissueconc, fitr1, fitdce, runmodel,
  defaultdict, ccc, nlsfit, demo, makeplots, validate

const verbose = true
const version = "v0.3"

include("util.jl")
include("fitting.jl")
include("models.jl")
include("science.jl")
include("plotting.jl")
include("demo.jl")
include("validate.jl")

end
