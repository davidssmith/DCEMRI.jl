using DCEMRI

isdir("results") || mkdir("results")
println("Running analysis of QIBA v6 (noisefree) data ...")
opts = defaultdict()
opts["datafile"] = "qiba6.mat"
opts["outfile"] = "results/results.mat"
opts["modelflags"] = 2
results = runmodel(opts)

println("Plotting results ...")
include("plotresults6.jl")
makeplots6(results, "results", dx=10)

isdir("results_noisy") || mkdir("results_noisy")
println("Running analysis of QIBA v6 (noisy) data ...")
opts = defaultdict()
opts["datafile"] = "qiba6noisy.mat"
opts["outfile"] = "results_noisy/results.mat"
opts["modelflags"] = 2
results = runmodel(opts)

println("Plotting results ...")
include("plotresults6.jl")
makeplots6(results, "results_noisy")
