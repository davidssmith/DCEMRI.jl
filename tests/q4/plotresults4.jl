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

  # AIF
  clf()
  plot(mat["t"], mat["aif"], "ko-")
  xlabel("time (min)")
  yticks([0:5])
  ylabel("[Gd-DTPA] (mM)")
  title("arterial input function, \$C_p\$")
  savefig("$outdir/aif.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(modelmap, interpolation="nearest", cmap="cubehelix")
  title("model used")
  colorbar(ticks=[0,1,2,3])
  savefig("$outdir/modelmap.pdf")

  # PARAMETER MAPS
  figure(figsize=(3, 6))
  clf()
  imshow(Kt, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.2)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  colorbar()
  savefig("$outdir/Kt.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_e\$")
  colorbar()
  savefig("$outdir/ve.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(vp, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.12)
  title("\$v_p\$")
  colorbar()
  savefig("$outdir/vp.pdf")

  figure(figsize=(3, 6))
  clf()
  imshow(resid, interpolation="nearest", cmap="cubehelix", vmin=0)
  title("residual")
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
  colorbar()
  savefig("$outdir/Kt_error.pdf")

  figure(figsize=(3, 6))
  clf()
  m = maximum(abs(ve_error))
  imshow(ve_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_e\$")
  colorbar()
  savefig("$outdir/ve_error.pdf")

  figure(figsize=(3, 6))
  clf()
  m = maximum(abs(vp_error))
  imshow(vp_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_p\$")
  colorbar()
  savefig("$outdir/vp_error.pdf")

end
