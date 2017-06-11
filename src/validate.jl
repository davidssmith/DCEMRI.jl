using PyPlot

function analyzer(mat::Dict, outdir::AbstractString; dx=1, makeplots=true, isExt=false)

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

  # compare to known truths
  if ~isExt
    u = ones(div(10,dx),div(10,dx))
    Kt_truth = repmat([0.01*u; 0.02*u; 0.05*u; 0.1*u; 0.2*u; 0.35*u], 1, 5)
    ve_truth = repmat([0.01*u; 0.05*u; 0.1*u; 0.2*u; 0.5*u]', 6, 1)
  else
    u = ones(div(180,dx),div(10,dx))
    Kt_truth = hcat(0.01u, 0.02u, 0.05u, 0.1u, 0.2u)

    u = ones(div(10,dx),div(50,dx))
    u = vcat(0.1u, 0.2u, 0.5u)
    ve_truth = vcat(u, u, u, u, u, u)

    u = ones(div(30,dx),div(50,dx))
    vp_truth = vcat(0.001u, 0.005u, 0.01u, 0.02u, 0.05u, 0.1u)
  end

  Kt_error = clamp.(100.0*(Kt - Kt_truth) ./ (Kt_truth + eps()), -100.0, 100.0)
  ve_error = clamp.(100.0*(ve - ve_truth) ./ (ve_truth + eps()), -100.0, 100.0)
  cccKt = ccc(Kt_truth, Kt)
  cccve = ccc(ve_truth, ve)
  cccvp = NaN
  print_with_color(:green, "Kt\n\tRMSE:\t$(sqrt(norm(Kt_error)^2 / length(Kt_error))) %\n")
  print_with_color(:green, "\terrmax:\t$(maximum(abs.(Kt_error)))\n")
  print_with_color(:green, "\tCCC:\t$cccKt\n")
  print_with_color(:green, "ve\n\tRMSE:\t$(sqrt(norm(ve_error)^2 / length(ve_error))) %\n")
  print_with_color(:green, "\terrmax:\t$(maximum(abs.(ve_error)))\n")
  print_with_color(:green, "\tCCC:\t$cccve\n")

  if isExt
    vp_error = clamp.(100.0*(vp - vp_truth) ./ (vp_truth + eps()), -100.0, 100.0)
    cccvp = ccc(vp_truth, vp)
    print_with_color(:green, "vp\n\tRMSE:\t$(sqrt(norm(vp_error)^2 / length(vp_error))) %\n")
    print_with_color(:green, "\terrmax:\t$(maximum(abs.(vp_error))) %\n")
    print_with_color(:green, "\tCCC:\t$cccvp\n")
  end

  if !makeplots
    return (cccKt, cccve, cccvp)
  end

  if ~isExt
    ytpos = collect((0+floor(Integer, 5/dx)):div(10,dx):(div(60,dx)-1))
    xtpos = collect((0+floor(Integer, 5/dx)):div(10,dx):(div(50,dx)-1))
    ytlabels = [string(x) for x in [0.01,0.02,0.05,0.1,0.2,0.35]]
    xtlabels = [string(x) for x in [0.01,0.05,0.1,0.2,0.5]]
    mapWidth = 4.5
    mapHeight = 4.5
    mapLabelX = "\$v_\\mathrm{e}\$"
    mapLabelY = "\$K^\\mathrm{trans}\$"
  else
    ytpos = collect((div(10,dx)+floor(Integer, 5/dx)):div(30,dx):(div(180,dx)-1))
    xtpos = collect((0+floor(Integer, 5/dx)):div(10,dx):(div(50,dx)-1))
    ytlabels = [string(x) for x in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]]
    xtlabels = [string(x) for x in [0.01,0.02,0.05,0.1,0.2]]
    mapWidth = 3.5
    mapHeight = 6
    mapLabelX = "\$K^\\mathrm{trans}\$"
    mapLabelY = "\$v_\\mathrm{p}\$"
  end

  println("Plotting results ...")
  # AIF
  figure(figsize=(4.5,4.5))
  clf()
  plot(mat["t"], mat["Cp"], "ko-")
  xlabel("time (min)")
  yticks(collect(0:2:10)) # This produces an error
  ylim(0,10)
  ylabel("[Gd-DTPA] (mM)")
  title("arterial input function, \$C_p\$")
  savefig("$outdir/aif.pdf",bbox_inches="tight")

  figure(figsize=(mapWidth, mapHeight))
  clf()
  imshow(modelmap, interpolation="nearest", cmap="cubehelix")
  title("model used")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel(mapLabelX)
  ylabel(mapLabelY)
  colorbar()
  savefig("$outdir/modelmap.pdf",bbox_inches="tight")

  # PARAMETER MAPS
  figure(figsize=(mapWidth, mapHeight))
  clf()
  imshow(Kt, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.35)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel(mapLabelX)
  ylabel(mapLabelY)
  colorbar()
  savefig("$outdir/Kt.pdf",bbox_inches="tight")

  figure(figsize=(mapWidth, mapHeight))
  clf()
  imshow(ve, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.6)
  title("\$v_\\mathrm{e}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel(mapLabelX)
  ylabel(mapLabelY)
  colorbar()
  savefig("$outdir/ve.pdf",bbox_inches="tight")

  if isExt
    figure(figsize=(mapWidth, mapHeight))
    clf()
    imshow(vp, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=0.12)
    title("\$v_p\$")
    xticks(xtpos, xtlabels, fontsize=8)
    yticks(ytpos, ytlabels)
    xlabel(mapLabelX)
    ylabel(mapLabelY)
    colorbar()
    savefig("$outdir/vp.pdf",bbox_inches="tight")
  end

  figure(figsize=(mapWidth, mapHeight))
  clf()
  imshow(resid, interpolation="nearest", cmap="cubehelix", vmin=0)
  title("residual")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel(mapLabelX)
  ylabel(mapLabelY)
  colorbar()
  savefig("$outdir/resid.pdf",bbox_inches="tight")

  figure(figsize=(mapWidth, mapHeight))
  clf()
  m = maximum(abs.(Kt_error))
  imshow(Kt_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$K^\\mathrm{trans}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel(mapLabelX)
  ylabel(mapLabelY)
  colorbar()
  savefig("$outdir/Kt_error.pdf",bbox_inches="tight")

  figure(figsize=(mapWidth, mapHeight))
  clf()
  m = maximum(abs.(ve_error))
  imshow(ve_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
  title("% error in \$v_\\mathrm{e}\$")
  xticks(xtpos, xtlabels)
  yticks(ytpos, ytlabels)
  xlabel(mapLabelX)
  ylabel(mapLabelY)
  colorbar()
  savefig("$outdir/ve_error.pdf",bbox_inches="tight")

  if isExt
    figure(figsize=(mapWidth, mapHeight))
    clf()
    m = maximum(abs.(vp_error))
    imshow(vp_error, interpolation="nearest", cmap="PiYG", vmin=-m, vmax=m)
    title("% error in \$v_p\$")
    xticks(xtpos, xtlabels, fontsize=8)
    yticks(ytpos, ytlabels)
    xlabel(mapLabelX)
    ylabel(mapLabelY)
    colorbar()
    savefig("$outdir/vp_error.pdf",bbox_inches="tight")
  end

  return (cccKt, cccve, cccvp)
end

function analyze(n, mat::Dict, outdir::AbstractString; kwargs...)
  if n == 4
    return analyzer(mat, outdir; isExt=true, kwargs...)
  elseif n == 6
    return analyzer(mat, outdir; isExt=false, kwargs...)
  end
end


function validate(n, outdir::AbstractString; kwargs...)
  @assert n == 4 || n == 6 "n must be 4 or 6"
  cd(Pkg.dir("DCEMRI/test/q$n"))

  println("Running analysis of noise-free QIBA v$n data ...")
  isdir("$outdir/results") || mkdir("$outdir/results")
  results = fitdata(datafile="qiba$n.mat",outfile="$outdir/results/results.mat",save=false)
  ccc = analyze(n, results, "$outdir/results", dx=10; kwargs...)

  println("Running analysis of noisy QIBA v$n data ...")
  isdir("$outdir/results_noisy") || mkdir("$outdir/results_noisy")
  results = fitdata(datafile="qiba$(n)noisy.mat", outfile="$outdir/results_noisy/results.mat",save=false)

  cccnoisy = analyze(n, results, "$outdir/results_noisy"; kwargs...)
  println("Validation complete. If saved, results can be found in $outdir.")
  return (ccc, cccnoisy)
end

validate(n; kwargs...) = validate(n, Pkg.dir("DCEMRI/test/q$n"); kwargs...)
function validate(kwargs...)
  ccc6, cccnoisy6 = validate(6; kwargs...)
  ccc4, cccnoisy4 = validate(4; kwargs...)
end
