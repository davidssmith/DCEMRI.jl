using PyPlot

function makeplots6(mat::Dict, outdir::AbstractString; dx=1)

  R1map = mat["R10"]
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

  ytpos = collect((0+floor(Integer, 5/dx)):div(10,dx):(div(60,dx)-1))
  xtpos = collect((0+floor(Integer, 5/dx)):div(10,dx):(div(50,dx)-1))
  ytlabels = [string(x) for x in [0.01,0.02,0.05,0.1,0.2,0.35]]
  xtlabels = [string(x) for x in [0.01,0.05,0.1,0.2,0.5]]

  # AIF
  figure(figsize=(4.5,4.5))
  clf()
  plot(mat["t"], mat["Cp"], "ko-")
  xlabel("time (min)")
  # yticks([0:2:10]) # This produces an error
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
  savefig("$outdir/modelmap.pdf",bbox_inches="tight")

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
  savefig("$outdir/Kt.pdf",bbox_inches="tight")

  figure(figsize=(4.5, 4))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_\\mathrm{e}\$")
  colorbar(ticks=[0:6]/10.0)
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  savefig("$outdir/ve.pdf",bbox_inches="tight")

  figure(figsize=(4.5, 4))
  clf()
  imshow(resid, interpolation="nearest", cmap="cubehelix", vmin=0)
  title("residual")
  colorbar()
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  savefig("$outdir/resid.pdf",bbox_inches="tight")

  # compare to known truths
  u = ones(div(10,dx),div(10,dx))
  Kt_truth = repmat([0.01*u; 0.02*u; 0.05*u; 0.1*u; 0.2*u; 0.35*u], 1, 5)
  ve_truth = repmat([0.01*u; 0.05*u; 0.1*u; 0.2*u; 0.5*u]', 6, 1)

  Kt_error = clamp.(100.0*(Kt - Kt_truth) ./ (Kt_truth + eps()), -100.0, 100.0)
  ve_error = clamp.(100.0*(ve - ve_truth) ./ (ve_truth + eps()), -100.0, 100.0)
  println("Kt RMSE (%): ", sqrt(norm(Kt_error)^2 / length(Kt_error)))
  println("Kt max error: ", maximum(abs.(Kt_error)))
  println("Kt CCC: ", ccc(Kt_truth, Kt))
  println("ve RMSE (%): ", sqrt(norm(ve_error)^2 / length(ve_error)))
  println("ve max error: ", maximum(abs.(ve_error)))
  println("ve CCC: ", ccc(ve_truth, ve))

  figure(figsize=(4.5, 4))
  clf()
  m = maximum(abs.(Kt_error))
  imshow(Kt_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$K^\\mathrm{trans}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar()
  savefig("$outdir/Kt_error.pdf",bbox_inches="tight")

  figure(figsize=(4.5, 4))
  clf()
  m = maximum(abs.(ve_error))
  imshow(ve_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_\\mathrm{e}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar()
  savefig("$outdir/ve_error.pdf",bbox_inches="tight")
end



function makeplots4(mat::Dict, outdir::AbstractString; dx=1)
  R1map = mat["R10"]
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

  ytpos = collect((div(10,dx)+floor(Integer, 5/dx)):div(30,dx):(div(180,dx)-1))
  xtpos = collect((0+floor(Integer, 5/dx)):div(10,dx):(div(50,dx)-1))
  ytlabels = [string(x) for x in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]]
  xtlabels = [string(x) for x in [0.01,0.02,0.05,0.1,0.2]]

  # AIF
  figure(figsize=(4,4))
  clf()
  plot(mat["t"], mat["Cp"], "ko-")
  xlabel("time (min)")
  #yticks([0:10])
  ylabel("[Gd-DTPA] (mM)")
  title("arterial input function, \$C_p\$")
  savefig("$outdir/aif.pdf",bbox_inches="tight")

  figure(figsize=(3.5, 6))
  clf()
  imshow(modelmap, interpolation="nearest", cmap="cubehelix")
  title("model used")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar(ticks=[0,1,2,3])
  savefig("$outdir/modelmap.pdf",bbox_inches="tight")

  # PARAMETER MAPS
  figure(figsize=(3.5, 6))
  clf()
  imshow(Kt, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.2)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/Kt.pdf",bbox_inches="tight")

  figure(figsize=(3.5, 6))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_e\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/ve.pdf",bbox_inches="tight")

  figure(figsize=(3.5, 6))
  clf()
  imshow(vp, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.12)
  title("\$v_p\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/vp.pdf",bbox_inches="tight")

  figure(figsize=(3.5, 6))
  clf()
  imshow(resid, interpolation="nearest", cmap="cubehelix", vmin=0)
  title("residual")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/resid.pdf",bbox_inches="tight")

  # compare to known truths
  u = ones(div(180,dx),div(10,dx))
  Kt_truth = hcat(0.01u, 0.02u, 0.05u, 0.1u, 0.2u)

  u = ones(div(10,dx),div(50,dx))
  u = vcat(0.1u, 0.2u, 0.5u)
  ve_truth = vcat(u, u, u, u, u, u)

  u = ones(div(30,dx),div(50,dx))
  vp_truth = vcat(0.001u, 0.005u, 0.01u, 0.02u, 0.05u, 0.1u)

  Kt_error = clamp.(100.0*(Kt - Kt_truth) ./ (Kt_truth + eps()), -100.0, 100.0)
  ve_error = clamp.(100.0*(ve - ve_truth) ./ (ve_truth + eps()), -100.0, 100.0)
  vp_error = clamp.(100.0*(vp - vp_truth) ./ (vp_truth + eps()), -100.0, 100.0)
  println("Kt RMSE (%): ", sqrt(norm(Kt_error)^2 / length(Kt_error)))
  println("Kt max error (%): ", maximum(abs.(Kt_error)))
  println("Kt CCC: ", ccc(Kt_truth, Kt))
  println("ve RMSE (%): ", sqrt(norm(ve_error)^2 / length(ve_error)))
  println("ve max error (%): ", maximum(abs.(ve_error)))
  println("ve CCC: ", ccc(ve_truth, ve))
  println("vp RMSE (%): ", sqrt(norm(vp_error)^2 / length(vp_error)))
  println("vp max error (%): ", maximum(abs.(vp_error)))
  println("vp CCC: ", ccc(vp_truth, vp))

  figure(figsize=(3.5, 6))
  clf()
  m = maximum(abs.(Kt_error))
  imshow(Kt_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$K^\\mathrm{trans}\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/Kt_error.pdf",bbox_inches="tight")

  figure(figsize=(3.5, 6))
  clf()
  m = maximum(abs.(ve_error))
  imshow(ve_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_e\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/ve_error.pdf",bbox_inches="tight")

  figure(figsize=(3.5, 6))
  clf()
  m = maximum(abs.(vp_error))
  imshow(vp_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_p\$")
  xticks(xtpos, xtlabels, fontsize=8)
  yticks(ytpos, ytlabels)
  xlabel("\$K^\\mathrm{trans}\$")
  ylabel("\$v_\\mathrm{p}\$")
  colorbar()
  savefig("$outdir/vp_error.pdf",bbox_inches="tight")
end

function makeplots(n, mat::Dict, outdir::AbstractString; dx=1)
  if n == 4
    makeplots4(mat, outdir; dx=dx)
  elseif n == 6
    makeplots6(mat, outdir; dx=dx)
  end
end


function validate(n, outdir::AbstractString)
  @assert n == 4 || n == 6 "n must be 4 or 6"
  cd(Pkg.dir("DCEMRI/test/q$n"))

  println("Running analysis of noise-free QIBA v$n data ...")
  isdir("$outdir/results") || mkdir("$outdir/results")
  results = fitdata(datafile="qiba$n.mat",outfile="$outdir/results/results.mat")
  println("Plotting results ...")
  makeplots(n, results, "$outdir/results", dx=10)

  println("Running analysis of noisy QIBA v$n data ...")
  isdir("$outdir/results_noisy") || mkdir("$outdir/results_noisy")
  results = fitdata(datafile="qiba$(n)noisy.mat",
                     outfile="$outdir/results_noisy/results.mat")

  println("Plotting results ...")
  makeplots(n, results, "$outdir/results_noisy")
  println("Validation complete. Results can be found in $outdir.")
end

validate(n) = validate(n, Pkg.dir("DCEMRI/test/q$n"))
function validate()
  validate(6)
  validate(4)
end
