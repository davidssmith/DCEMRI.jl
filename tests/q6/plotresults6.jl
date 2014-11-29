using PyPlot
using DCEMRI

function makeplots(mat::Dict, outdir::String; dx::Int64=1)

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

  figure(figsize=(4.5, 4))
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
  figure(figsize=(4.5, 4))
  clf()
  imshow(Kt, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.35)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  colorbar(ticks=[0.1,0.2,0.3,0.4])
  savefig("$outdir/Kt.pdf")

  figure(figsize=(4.5, 4))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_\\mathrm{e}\$")
  colorbar(ticks=[0:6]/10.0)
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel("\$v_\\mathrm{e}\$")
  ylabel("\$K^\\mathrm{trans}\$")
  savefig("$outdir/ve.pdf")

  figure(figsize=(4.5, 4))
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

  figure(figsize=(4.5, 4))
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

  figure(figsize=(4.5, 4))
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

#  clf()
#  #x = [[x,y] for x in Kt_truth[:], y in  Kt[:]]
#  x = mean(Kt_truth, 2)
#  y = mean(Kt, 2)
#  yerr = 2*std(Kt, 2)
#  plot([0,0.4],[0,0.4], "k:")
#  #plot(Kt_truth[:] + 0.003*randn(length(Kt_truth[:])), Kt[:], "k.", alpha=0.2)
#  errorbar(x, y, yerr=yerr, fmt="ko")
#  xlabel("true \$K^\\mathrm{trans}\$")
#  ylabel("derived \$K^\\mathrm{trans}\$")
#  xlim(0,0.4)
#  ylim(0,0.4)
#  savefig("$outdir/kt_boxplot.pdf")

#  clf()
#  x = mean(ve_truth, 1).'
#  y = mean(ve, 1).'
#  yerr = 2*std(ve, 1).'
#  plot([0,1],[0,1], "k:")
#  #plot(ve_truth[:] + 0.003*randn(length(ve_truth[:])), ve[:],"k.", alpha=0.2)
#  errorbar(x, y, yerr=yerr, fmt="ko")
#  xlabel("true \$v_\\mathrm{e}\$")
#  ylabel("derived \$v_\\mathrm{e}\$")
#  xlim(0,1)
#  xlim(0,1)
#  savefig("$outdir/ve_boxplot.pdf")

end
