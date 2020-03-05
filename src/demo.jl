function demo(outdir::AbstractString="results")
  cd(joinpath(dirname(pathof(DCEMRI)), "..", "demo"))
  isdir(outdir) || mkdir(outdir)
  println("Processing in vivo data ...")

  # run the model
  results = fitdata(datafile="invivo.mat", outfile="$outdir/results.mat", models=[2])

  if !haskey(Pkg.installed(), "PyPlot")
    # Do no make plots if PyPlot not installed
    println("PyPlot not installed. Plots will not be produced.")
  else
    # plot the results
    println("Plotting results ...")
    makeplots(results; outdir=outdir)
    if outdir == "results"
      println("Results can be found in ", Pkg.dir("DCEMRI/demo/$outdir"))
    else
      println("Results can be found in $outdir")
    end
  end

  println("Demo run complete.")
end
