cd("q4")
if !isfile("qiba4.mat")
  run(`python prep4.py`)
end
include("q4/validate4.jl")

cd("../q6")
if !isfile("qiba6.mat")
  run(`python prep6.py`)
end
include("q6/validate6.jl")

cd("..")
