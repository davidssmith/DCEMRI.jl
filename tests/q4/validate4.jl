using DCEMRI

println("Running analysis of QIBA v4 extended Tofts data ...")
opts = defaultdict()
opts["datafile"] = "qiba4.mat"
isdir("results") || mkdir("results")
opts["outfile"] = "results/results.mat"
opts["modelflags"] = 4
delete!(opts, "mask")
results = runmodel(opts)

println("Plotting results ...")
include("plotresults4.jl")
makeplots(results, "results", dx=10)


println("Running analysis of noisy QIBA v4 extended Tofts data ...")
opts = defaultdict()
opts["datafile"] = "qiba4noisy.mat"
isdir("results_noisy") || mkdir("results_noisy")
opts["outfile"] = "results_noisy/results.mat"
opts["modelflags"] = 4
delete!(opts, "mask")
results = runmodel(opts)

println("Plotting results ...")
makeplots(results, "results_noisy")
