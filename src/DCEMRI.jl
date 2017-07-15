module DCEMRI

using ArgParse
using Calculus    #.jacobian
using MAT

export ser, r1eff, tissueconc, fitr1, fitdce, fitdata,
  defaults, ccc, nlsfit, makeplots, demo, validate

const verbose = true
const version = v"0.1.2"

if Pkg.installed("PyPlot")==Void()
  # println("Optional package (PyPlot) not installed.")
else
  using PyPlot
end

include("util.jl")
include("fitting.jl")
include("models.jl")
include("science.jl")
include("plotting.jl")
include("demo.jl")
include("validate.jl")

end
