

function makeplots6(mat::Dict, outdir::String; dx::Int64=1)

  R1map = mat["R1"]
  S0map = mat["S0"]
  modelmap = mat["modelmap"]
  Ct = mat["Ct"]
  Kt = mat["Kt"]
  ve = mat["ve"]
  vp = mat["vp"]
  resid = mat["resid"]
  q = quantile(S0map[:], 0.99)
  S0map[S0map .> q] = q
  back = (S0map - minimum(S0map)) / (maximum(S0map) - minimum(S0map))
  mask = convert(Array{Bool,2}, mat["mask"])

  ytpos = [(0+ifloor(5/dx)):div(10,dx):(div(60,dx)-1)]
  xtpos = [(0+ifloor(5/dx)):div(10,dx):(div(50,dx)-1)]
  ytlabels = [string(x) for x in [0.01,0.02,0.05,0.1,0.2,0.35]]
  xtlabels = [string(x) for x in [0.01,0.05,0.1,0.2,0.5]]

  # AIF
  figure(figsize=(4.5,4.5))
  clf()
  plot(mat["t"], mat["aif"], "ko-")
  xlabel("time (min)")
  yticks([0:2:10])
  ylim(0,10)
  ylabel("[Gd-DTPA] (mM)")
  title("arterial input function, \$C_p\$")
  savefig("$outdir/aif.pdf")

  figure(figsize=(4.5, 4.5))
  clf()
  imshow(modelmap, interpolation="nearest", cmap="cubehelix")
  title("model used")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar(ticks=[0,1,2,3])
  savefig("$outdir/modelmap.pdf")

  # PARAMETER MAPS
  figure(figsize=(4.5, 4.5))
  clf()
  imshow(Kt, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.35)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar(ticks=[0.1,0.2,0.3,0.4])
  savefig("$outdir/Kt.pdf")

  figure(figsize=(4.5, 4.5))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_\\mathrm{e}\$")
  colorbar(ticks=[0:6]/10.0)
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  savefig("$outdir/ve.pdf")

  figure(figsize=(4.5, 4.5))
  clf()
  imshow(resid, interpolation="nearest", cmap="cubehelix", vmin=0)
  title("residual")
  colorbar()
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  savefig("$outdir/resid.pdf")

  # compare to known truths
  u = ones(div(10,dx),div(10,dx))
  Kt_truth = repmat([0.01*u, 0.02*u, 0.05*u, 0.1*u, 0.2*u, 0.35*u], 1, 5)
  ve_truth = repmat([0.01*u, 0.05*u, 0.1*u, 0.2*u, 0.5*u]', 6, 1)

  Kt_error = clamp(100.0*(Kt - Kt_truth) ./ (Kt_truth + eps()), -100.0, 100.0)
  ve_error = clamp(100.0*(ve - ve_truth) ./ (ve_truth + eps()), -100.0, 100.0)
  println("Kt RMSE (%): ", sqrt(norm(Kt_error)^2 / length(Kt_error)))
  println("Kt max error: ", maximum(abs(Kt_error)))
  println("Kt CCC: ", ccc(Kt_truth, Kt))
  println("ve RMSE (%): ", sqrt(norm(ve_error)^2 / length(ve_error)))
  println("ve max error: ", maximum(abs(ve_error)))
  println("ve CCC: ", ccc(ve_truth, ve))

  figure(figsize=(4.5, 4.5))
  clf()
  m = maximum(abs(Kt_error))
  imshow(Kt_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$K^\\mathrm{trans}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar()
  savefig("$outdir/Kt_error.pdf")

  figure(figsize=(4.5, 4.5))
  clf()
  m = maximum(abs(ve_error))
  imshow(ve_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_\\mathrm{e}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar()
  savefig("$outdir/ve_error.pdf")
end



function makeplots4(mat::Dict, outdir::String; dx::Int64=1)
  R1map = mat["R1"]
  S0map = mat["S0"]
  modelmap = mat["modelmap"]
  Ct = mat["Ct"]
  Kt = mat["Kt"]
  ve = mat["ve"]
  vp = mat["vp"]
  resid = mat["resid"]
  q = quantile(S0map[:], 0.99)
  S0map[S0map .> q] = q
  back = (S0map - minimum(S0map)) / (maximum(S0map) - minimum(S0map))
  mask = convert(Array{Bool,2}, mat["mask"])

  ytpos = [(div(10,dx)+ifloor(5/dx)):div(30,dx):(div(180,dx)-1)]
  xtpos = [(0+ifloor(5/dx)):div(10,dx):(div(50,dx)-1)]
  ytlabels = [string(x) for x in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]]
  xtlabels = [string(x) for x in [0.01,0.02,0.05,0.1,0.2]]

  # AIF
  figure(figsize=(4,4))
  clf()
  plot(mat["t"], mat["aif"], "ko-")
  xlabel("time (min)")
  #yticks([0:10])
  ylabel("[Gd-DTPA] (mM)")
  title("arterial input function, \$C_p\$")
  savefig("$outdir/aif.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(modelmap, interpolation="nearest", cmap="cubehelix")
  title("model used")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar(ticks=[0,1,2,3])
  savefig("$outdir/modelmap.pdf")

  # PARAMETER MAPS
  figure(figsize=(3, 6))
  clf()
  imshow(Kt, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.2)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/Kt.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_e\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/ve.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(vp, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.12)
  title("\$v_p\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/vp.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(resid, interpolation="nearest", cmap="cubehelix", vmin=0)
  title("residual")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/resid.pdf")

  # compare to known truths
  u = ones(div(180,dx),div(10,dx))
  Kt_truth = hcat(0.01u, 0.02u, 0.05u, 0.1u, 0.2u)

  u = ones(div(10,dx),div(50,dx))
  u = vcat(0.1u, 0.2u, 0.5u)
  ve_truth = vcat(u, u, u, u, u, u)

  u = ones(div(30,dx),div(50,dx))
  vp_truth = vcat(0.001u, 0.005u, 0.01u, 0.02u, 0.05u, 0.1u)

  Kt_error = clamp(100.0*(Kt - Kt_truth) ./ (Kt_truth + eps()), -100.0, 100.0)
  ve_error = clamp(100.0*(ve - ve_truth) ./ (ve_truth + eps()), -100.0, 100.0)
  vp_error = clamp(100.0*(vp - vp_truth) ./ (vp_truth + eps()), -100.0, 100.0)
  println("Kt RMSE (%): ", sqrt(norm(Kt_error)^2 / length(Kt_error)))
  println("Kt max error (%): ", maximum(abs(Kt_error)))
  println("Kt CCC: ", ccc(Kt_truth, Kt))
  println("ve RMSE (%): ", sqrt(norm(ve_error)^2 / length(ve_error)))
  println("ve max error (%): ", maximum(abs(ve_error)))
  println("ve CCC: ", ccc(ve_truth, ve))
  println("vp RMSE (%): ", sqrt(norm(vp_error)^2 / length(vp_error)))
  println("vp max error (%): ", maximum(abs(vp_error)))
  println("vp CCC: ", ccc(vp_truth, vp))

  figure(figsize=(3, 6))
  clf()
  m = maximum(abs(Kt_error))
  imshow(Kt_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$K^\\mathrm{trans}\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/Kt_error.pdf")

  figure(figsize=(3, 6))
  clf()
  m = maximum(abs(ve_error))
  imshow(ve_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_e\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/ve_error.pdf")

  figure(figsize=(3, 6))
  clf()
  m = maximum(abs(vp_error))
  imshow(vp_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_p\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/vp_error.pdf")
end


function validate4(outdir::String)
  cd(Pkg.dir("DCEMRI/tests/q4"))
  println("Running analysis of QIBA v4 extended Tofts data ...")
  opts = defaultdict()
  opts["datafile"] = "qiba4.mat"
  isdir("$outdir/results") || mkdir("$outdir/results")
  opts["outfile"] = "$outdir/results/results.mat"
  opts["modelflags"] = 4
  delete!(opts, "mask")
  results = runmodel(opts)

  println("Plotting results ...")
  makeplots4(results, "$outdir/results", dx=10)

  println("Running analysis of noisy QIBA v4 extended Tofts data ...")
  opts = defaultdict()
  opts["datafile"] = "qiba4noisy.mat"
  isdir("$outdir/results_noisy") || mkdir("$outdir/results_noisy")
  opts["outfile"] = "$outdir/results_noisy/results.mat"
  opts["modelflags"] = 4
  delete!(opts, "mask")
  results = runmodel(opts)

  println("Plotting results ...")
  makeplots4(results, "$outdir/results_noisy")
end
validate4() = validate4(Pkg.dir("DCEMRI/tests/q4"))

function validate6(outdir::String)
  cd(Pkg.dir("DCEMRI/tests/q6"))
  println("Running analysis of QIBA v6 standard Tofts data ...")
  opts = defaultdict()
  opts["datafile"] = "qiba6.mat"
  isdir("$outdir/results") || mkdir("$outdir/results")
  opts["outfile"] = "$outdir/results/results.mat"
  opts["modelflags"] = 2
  delete!(opts, "mask")
  results = runmodel(opts)

  println("Plotting results ...")
  makeplots6(results, "$outdir/results", dx=10)

  println("Running analysis of noisy QIBA v6 standard Tofts data ...")
  opts = defaultdict()
  opts["datafile"] = "qiba6noisy.mat"
  isdir("$outdir/results_noisy") || mkdir("$outdir/results_noisy")
  opts["outfile"] = "$outdir/results_noisy/results.mat"
  opts["modelflags"] = 2
  delete!(opts, "mask")
  results = runmodel(opts)

  println("Plotting results ...")
  makeplots6(results, "$outdir/results_noisy")
end
validate6() = validate6(Pkg.dir("DCEMRI/tests/q6"))

function validate()
  validate6()
  validate4()
end
