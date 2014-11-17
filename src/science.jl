function ser{N}(images::Array{Float64,N}, thresh::Float64=0.01)
  @dprint "computing signal enhancement ratios"
  @assert thresh > 0.0
  dims = size(images)
  nt = dims[1]
  n = prod(dims[2:end])
  SER = zeros(n)
  images = reshape(images, (nt, n))
  thresh *= maximum(images[:])
  for k in 1:n
    SER[k] = sum(images[1:3,k]) > thresh ? sum(images[end-3:end,k]) / sum(images[1:3,k]) : 0.0
  end
  reshape(SER, dims[2:end])
end

function r1eff(S::Array{Float64,3}, S0::Matrix{Float64}, R1::Matrix{Float64},
               TR::Float64, flip::Float64)
  @dprint "converting DCE signal to effective R1"
  @assert 0.0 < flip "flip angle must be positive"
  @assert 0.0 < TR && TR < 1.0 "TR must be in units of ms"
  nt, nx, ny = size(S)
  A = copy(S)
  R1eff = similar(S)
  for x in 1:nx, y in 1:ny
    E0 = exp(-R1[x,y] * TR)
    A[:,x,y] = A[:,x,y] / S0[x,y] # normalize by pre-contrast signal
    for t in 1:nt
      E = (1.0 - A[t,x,y] + A[t,x,y]*E0 - E0*cos(flip)) /
        (1.0 - A[t,x,y]*cos(flip) + A[t,x,y]*E0.*cos(flip) - E0*cos(flip))
      R1eff[t,x,y] = E > 0.0 ? (-1.0 / TR) * log(E) : 0.0
    end
  end
  R1eff
end

function tissueconc{M,N}(R1eff::Array{Float64,M}, R1map::Array{Float64,N}, relaxivity::Float64)
  @dprint "converting effective R1 to tracer tissue concentration Ct"
  @assert relaxivity > 0.0
  Ct = similar(R1eff)
  nt = size(R1eff,1)
  xidxs = find(R1map)
  for x in xidxs, t in 1:nt
    Ct[t,x] = R1eff[t,x] > 0.0 ? (R1eff[t,x] - R1map[x]) / relaxivity : 0.0
  end
  Ct
end


function fitr1(images, flip_angles::Vector{Float64}, TR::Float64,
               resid_thresh::Float64=0.01)
  @dprint "fitting R1 relaxation rate to multi-flip data"
  sizein = size(images)
  n = prod(sizein[2:end])
  nangles = sizein[1]
  @assert nangles == length(flip_angles)
  images = reshape(images, (nangles, n))
  p0 = [maximum(images), 1.0]
  model(x,p) = spgreqn(x, p, TR)
  idxs = find(mean(images, 1) .> 0.1*maximum(images))
  params, resid = nlsfit(model, images, idxs, flip_angles, p0)
  S0map = reshape(params[1,:], sizein[2:end])
  R1map = reshape(params[2,:], sizein[2:end])
  (R1map, S0map, resid)
end


function fitdce{N}(Ct::Array{Float64,N}, mask::BitMatrix, t::Vector{Float64},
                   Cp::Vector{Float64}; modelflags::Int=7, verbose::Bool=false,
                   residthresh::Float64=1.0, ktcutoff::Float64=5.0)
  @dprint "fitting DCE data"
  sizein = size(Ct)
  n = prod(sizein[2:end])
  nt = sizein[1]
  @assert nt == length(t)
  Ct = reshape(Ct, (nt, n))
  idxs = find(mask)

  nidxs = length(idxs)
  nmodels = sum(int(split(bits(modelflags),"")))
  @assert nmodels > 0 "at least one model must be specified in modelflags"
  resid = Inf*ones(n)
  params = zeros(3, n)
  modelmap = zeros(Uint8, n)

  if modelflags & 4 != 0
    @dprint "attempting Extended Tofts-Kety model"
    p0 = [0.01, 0.01, 0.01]
    f3(x,p) = extendedtoftskety(x, p, Cp)
    p, r, dof = nlsfit(f3, Ct, idxs, t, p0)
    p[2,idxs] = p[1,idxs] ./ p[2,idxs]
    resid[:] = squeeze(sumabs2(r, 1), 1) / dof
    params[:] = p
    modelmap[idxs] = 3
  end
  if modelflags & 2 != 0
    @dprint "attempting Standard Tofts-Kety model"
    p0 = [0.01, 0.01]
    f2(x,p) = toftskety(x, p, Cp)
    p, r, dof = nlsfit(f2, Ct, idxs, t, p0)
    r = squeeze(sumabs2(r, 1), 1) / dof
    for k in idxs
      if r[k] <= resid[k]
        resid[k] = r[k]
        params[1,k] = p[1,k]
        params[2,k] = p[1,k] / p[2,k]  # ve = Ktrans / kep
        params[3,k] = 0.0
        modelmap[k] = 2
      end
    end
  end
  if modelflags & 1 != 0
    @dprint "attempting plasma-only model"
    p0 = [0.01]
    f1(x,p) = onecompartment(x, p, Cp)
    p, r, dof = nlsfit(f1, Ct, idxs, t, p0)
    r = squeeze(sumabs2(r, 1), 1) / dof
    for k in idxs
      if r[k] <= resid[k]
        resid[k] = r[k]
        params[3,k] = p[1,k]
        params[1:2,k] = 0.0
        modelmap[k] = 1
      end
    end
  end
  for k in idxs
    if resid[k] <= residthresh
      params[1,k] = clamp(params[1,k], eps(), 5.0)
      params[2:end,k] =  clamp(params[2:end,k], eps(), 1.0)
    else
      params[:,k] = 0.0
      mask[k] = false
    end
  end
  params = reshape(params, [size(params,1), sizein[2:end]...]...)
  resid = reshape(resid, sizein[2:end])
  modelmap = reshape(modelmap, sizein[2:end])
  (params, resid, modelmap)
end


function runmodel(opts::Dict)
  @dprint "reading input data"

  # load DCE and R1 data
  mat = matread(opts["datafile"])
  validate(mat)
  const relaxivity = opts["relaxivity"] == nothing ? mat["R"] : opts["relaxivity"]
  const TR = opts["TR"] == nothing ? mat["TR"] : opts["TR"]
  const modelflags = haskey(mat, "modelflags") ? int(mat["modelflags"]) : opts["modelflags"]
  const outfile = haskey(mat, "outfile") ? mat["outfile"] : opts["outfile"]
  aifdata = mat["aif"][:]  # colon makes sure that we get an Array{T,1}

  # parallel workers are better than multithreaded BLAS for this problem
  # run julia with the '-p <n>' flag to start with n workers
  blas_set_num_threads(1)

  if haskey(mat, "R1map") && haskey(mat, "S0map")
    @dprint "found existing R1 map"
    R1map = mat["R1map"]
    S0map = mat["S0map"]
  else
    @dprint "found multi-flip data"
    t1data = mat["t1data"]
    t1flip = (isempty(opts["t1flip"]) ? mat["t1flip"][:] : opts["t1flip"]) * pi / 180.0
    @assert length(t1flip) == size(t1data,1) # must have one flip angle per T1 image
    R1map, S0map = fitr1(t1data, t1flip, TR)
  end

  dcedata = mat["dcedata"]
  t = mat["t"][:] / 60.0  # convert time to min
  dceflip = (opts["dceflip"] == nothing ? mat["dceflip"] : opts["dceflip"]) * pi / 180.0

  nt, nx, ny = size(dcedata)
  @assert nx == size(R1map,1) "R1 map and DCE images have different numbers of rows"
  @assert ny == size(R1map,2) "R1 map and DCE images have different numbers of columns"

  # MAIN: run postprocessing steps
  SER = ser(dcedata)
  S0 = squeeze(mean(dcedata[1:2,:,:],1),1)
  R1eff = r1eff(dcedata, S0, R1map, TR, dceflip)
  Ct = tissueconc(R1eff, R1map, relaxivity)
  mask = haskey(mat, "mask") ? convert(BitArray{2}, mat["mask"]) : SER .> 2.0
  params, resid, modelmap = fitdce(Ct, mask, t, aifdata, modelflags=modelflags,
                                   verbose=opts["verbose"])

  @dprint "saving results to $outfile"
  results = Dict()
  results["t"] = t
  results["aif"] = aifdata
  results["R1"] = R1map
  results["S0"] = S0map
  results["SER"] = SER
  results["mask"] = convert(Array{Bool,2}, mask)
  results["modelflags"] = modelflags
  results["modelmap"] = modelmap
  results["R1eff"] = R1eff
  results["Ct"] = Ct
  results["Kt"] = squeeze(params[1,:,:], 1)
  results["ve"] = squeeze(params[2,:,:], 1)
  results["vp"] = squeeze(params[3,:,:], 1)
  results["resid"] = resid
  matwrite(outfile, results)

  results
end

runmodel( ;kwargs...) = runmodel(kwargs2dict(kwargs))
