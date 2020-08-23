function jetrgb(x::Float64)
  y = 4x
  r = min(y - 1.5, -y + 4.5)
  g = min(y - 0.5, -y + 3.5)
  b = min(y + 0.5, -y + 2.5)
  [clamp(r,0.,1.), clamp(g,0.,1.), clamp(b,0.,1.)]
end

function grayrgb(x::Float64)
  y = clamp(x,0.,1.)
  [y,y,y]
end

jetgray(x::Float64) = x <= 0.01 ? grayrgb(100x) : jetrgb(x)

function redblue(x::Float64)
  cmap = colormap("RdBu")
  idx = round(Int,99x + 1)
  cmap[idx]
end

function cubehelix(x::Float64; start=0.3, rot=-0.5, gamma=1.0,
                   minSat=1.2, maxSat=1.2)
  angle = 2pi*(start/3.0 + rot*x + 1.)
  x = x.^gamma
  satar = x.*(maxSat - minSat) + minSat
  amp = satar.*x.*(1. - x) / 2.0
  red = x + amp.*(-0.14861*cos(angle) + 1.78277*sin(angle))
  grn = x + amp.*(-0.29227*cos(angle) - 0.90649*sin(angle))
  blu = x + amp.*(1.97294*cos(angle))
  red = clamp(red, 0., 1.)
  grn = clamp(grn, 0., 1.)
  blu = clamp(blu, 0., 1.)
  [red, grn, blu]
end

function overlay(front::Array{Float64,2}, back::Array{Float64,2}, mask::Array{Bool,2})
  @assert(length(front) == length(back))
  z = similar(front)
  scale= maximum(front)
  back[:] = scale*(back[:] - minimum(back)) / (maximum(back) - minimum(back))
  for k in 1:length(mask)
    z[k] = mask[k] ? max(front[k], nextfloat(nextfloat(0.01scale))) : 0.01back[k]
  end
  z
end

function oplot2(front::Array{Float64,2}, back::Array{Float64,2}, mask::Array{Bool,2})
  @assert length(front) == length(back)
  # normalize the images
  r = maximum(back) - minimum(back)
  if r == 0.0
    back[:] = 0.0
  else
    @. back = (back - minimum(back)) / r
  end
  n,m = size(front)
  img = zeros(n,m,3)
  for k in 1:m, j in 1:n
    if mask[j,k]
      s = front[j,k]
      s = clamp(s, 0.0, 1.0)
      cmidx = 100 - round(Int,99.0*s)
      img[j,k,:] = jetrgb(s)
    else
      s = back[j,k]
      if isnan(s)
        s = 0
      end
      cmidx = 100 - round(Int,99.0*s)
      img[j,k,1] = s
      img[j,k,2] = s
      img[j,k,3] = s
    end
  end
  img
end

function makeplots(mat::Dict; outdir::AbstractString="results")
  isdir(outdir) || mkdir(outdir)
  R1map = mat["R10"]
  S0map = mat["S0"]
  modelmap = mat["modelmap"]
  Ct = mat["Ct"]
  Kt = mat["Kt"]
  ve = mat["ve"]
  vp = mat["vp"]
  resid = mat["resid"]
  q = quantile(S0map[:], 0.99)
  S0map[S0map .> q] .= q
  back = @. (S0map - minimum(S0map)) / (maximum(S0map) - minimum(S0map))
  mask = convert(Array{Bool,2}, mat["mask"])

  figure(figsize=(4.5,4.5))
  clf()
  plot(mat["t"], mat["Cp"], "ko-")
  xlabel("time (min)")
  yticks(collect(0:5))
  ylabel("[Gd-DTPA] (mM)")
  title("arterial input function, \$C_p\$")
  savefig("$outdir/aif.pdf")

  figure()
  clf()
  imshow(mat["SER"], interpolation="nearest", cmap="cubehelix", vmin=0, vmax=10)
  colorbar()
  title("signal enhancement ratio")
  savefig("$outdir/ser.pdf")

  figure()
  clf()
  imshow(mask, interpolation="nearest", cmap="gray")
  title("mask")
  savefig("$outdir/mask.pdf")

  figure()
  clf()
  imshow(R1map, interpolation="nearest", cmap="cubehelix", vmin=0, vmax=5)
  title("\$R_1\$ relaxation rate (s\$^{-1}\$)")
  colorbar(ticks=collect(0:5))
  savefig("$outdir/R1.pdf")

  figure()
  clf()
  imshow(S0map, interpolation="nearest", cmap="gray")
  title("\$S_0\$")
  colorbar()
  savefig("$outdir/S0.pdf")

  figure()
  clf()
  Ct = dropdims(maximum(Ct; dims=1); dims=1)
  x = oplot2(clamp.(Ct, 0.0, 5.0), back, mask)
  imshow(x, interpolation="nearest", cmap="jet", vmin=0, vmax=5)
  title("max tissue conc., \$C_t\$ (mmol)")
  colorbar()
  savefig("$outdir/Ct.pdf")

  figure()
  clf()
  x = oplot2(clamp.(Kt, 0.0, 1.0), back, mask)
  imshow(x, interpolation="nearest", cmap="jet", vmin=0, vmax=1)
  title("\$K^\\mathrm{trans}\$ (min\$^{-1}\$)")
  colorbar(ticks=collect(0:2:10)/10.0)
  savefig("$outdir/Kt.pdf")

  figure()
  clf()
  x = oplot2(clamp.(ve, 0.0, 1.0), back, mask)
  imshow(x, interpolation="nearest", cmap="jet", vmin=0, vmax=1)
  title("\$v_e\$")
  colorbar(ticks=[0,0.2,0.4,0.6,0.8,1])
  savefig("$outdir/ve.pdf")

  figure()
  clf()
  x = oplot2(clamp.(vp, 0.0, 1.0), back, mask)
  imshow(x, interpolation="nearest", cmap="jet", vmin=0, vmax=1)
  title("\$v_p\$")
  colorbar(ticks=[0,0.2,0.4,0.6,0.8,1])
  savefig("$outdir/vp.pdf")

  figure()
  clf()
  kep = clamp.(Kt./ve, 0.0, 10.0)
  x = oplot2(kep, back, mask)
  imshow(x, interpolation="nearest", cmap="jet", vmin=0, vmax=10)
  title("\$k_{ep}\$")
  colorbar(ticks=[0,2,4,6,8,10])
  savefig("$outdir/kep.pdf")

  figure()
  clf()
  imshow(modelmap, interpolation="nearest", cmap="cubehelix")
  title("model used")
  colorbar(ticks=[0,1,2,3])
  savefig("$outdir/modelmap.pdf")

  figure()
  clf()
  x = oplot2(clamp.(100*resid, 0.0, 1.0), back, mask)
  imshow(x, interpolation="nearest", cmap="jet", vmin=0, vmax=1.0)
  title("residual \$\\times\$ 100")
  #colorbar(ticks=[0:2:10]/10000.0)
  colorbar()
  savefig("$outdir/resid.pdf")
end

function makeplots(matfile::AbstractString)
  println("creating plots from the file ", matfile)
  mat = read(matopen(matfile))
  makeplots(mat)
end
