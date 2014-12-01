function demo(outdir::String="results")
  cd(Pkg.dir("DCEMRI/demo"))
  isdir(outdir) || mkdir(outdir)
  println("Processing in vivo data ...")

  # run the model
  results = runmodel(datafile="invivo.mat", outfile="$outdir/results.mat", models=[2])

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
