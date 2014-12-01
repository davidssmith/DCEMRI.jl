function rundemo(outdir::String="results")
  cd(Pkg.dir("DCEMRI/demo"))
  isdir(outdir) || mkdir(outdir)
  println("Processing in vivo data ...")

  # initialize the model parameters with the defaults
  opts = defaultdict()

  # then modify as necessary to fit your data
  opts["datafile"] = "invivo.mat"
  opts["outfile"] = "$outdir/results.mat"
  opts["modelflags"] = 2

  # run the model
  results = runmodel(opts)

  # plot the results
  println("Plotting results ...")
  makeplots(results; outdir=outdir)

  println("Demo run complete.")
  if outdir == "results"
    println("Results can be found in ", Pkg.dir("DCEMRI/demo/$outdir"))
  else
    println("Results can be found in $outdir")
  end
end
