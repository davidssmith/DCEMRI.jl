# This in vivo example shows an example of what an interactive invocation of
# the module might look like.

using DCEMRI

isdir("results") || mkdir("results")
println("Processing in vivo data ...")

# initialize the model parameters with the defaults
opts = defaultdict()

# then modify as necessary to fit your data
opts["datafile"] = "invivo.mat"
opts["outfile"] = "results/results.mat"
opts["modelflags"] = 2

# run the model
results = runmodel(opts)

# plot the results
println("Plotting results ...")
include("../src/plotresults.jl")
makeplots(results)
