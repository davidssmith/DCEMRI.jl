function ser{N}(x::Array{Float64,N}, thresh::Float64=0.01)
  @dprint "computing signal enhancement ratios"
  @assert thresh > 0.0
  dims = size(x)
  nt = dims[1]
  n = prod(dims[2:end])
  SER = zeros(n)
  x = reshape(x, (nt, n))
  thresh *= maximum(x[:])
  for k in 1:n
    SER[k] = sum(x[1:3,k]) > thresh ? sum(x[end-3:end,k]) / sum(x[1:3,k]) : 0.0
  end
  reshape(SER, dims[2:end])
end

function r1eff{M,N}(S::Array{Float64,M}, R10::Array{Float64,N}, TR::Float64, flip::Float64)
  @dprint "converting DCE signal to effective R1"
  @assert 0.0 < flip "flip angle must be positive"
  @assert 0.0 < TR && TR < 1.0 "TR must be in units of ms"
  dims = size(S)
  nt = dims[1]
  n = prod(dims[2:end])
  S = reshape(S, (nt, n))
  S0 = mean(S[1:2,:],1)
  A = copy(S)
  R1 = similar(S)
  for k = 1:n
    E0 = exp(-R10[k] * TR)
    A[:,k] = A[:,k] / S0[k] # normalize by pre-contrast signal
    for t in 1:nt
      E = (1.0 - A[t,k] + A[t,k]*E0 - E0*cos(flip)) /
        (1.0 - A[t,k]*cos(flip) + A[t,k]*E0.*cos(flip) - E0*cos(flip))
      R1[t,k] = E > 0.0 ? (-1.0 / TR) * log(E) : 0.0
    end
  end
  reshape(R1, dims)
end

function tissueconc{M,N}(R1::Array{Float64,M}, R10::Array{Float64,N}, r1::Float64)
  @dprint "converting effective R1 to tracer tissue concentration Ct"
  @assert r1 > 0.0
  dims = size(R1)
  nt = dims[1]
  xidxs = find(R10)
  n = prod(dims[2:end])
  R1 = reshape(R1, (nt, n))
  Ct = similar(R1)
  for x in xidxs, t in 1:nt
    Ct[t,x] = R1[t,x] > 0.0 ? (R1[t,x] - R10[x]) / r1 : 0.0
  end
  R1 = reshape(R1, dims)
  reshape(Ct, dims)
end


function fitr1(x, flip_angles::Vector{Float64}, TR::Float64,
               resid_thresh::Float64=0.01)
  @dprint "fitting R1 relaxation rate to multi-flip data"
  sizein = size(x)
  n = prod(sizein[2:end])
  nangles = sizein[1]
  @assert nangles == length(flip_angles)
  x = reshape(x, (nangles, n))
  p0 = [maximum(x), 1.0]
  model(x,p) = spgreqn(x, p, TR)
  idxs = find(mean(x, 1) .> 0.1*maximum(x))
  params, resid = nlsfit(model, x, idxs, flip_angles, p0)
  S0 = reshape(params[1,:], sizein[2:end])
  R10 = reshape(params[2,:], sizein[2:end])
  (R10, S0, resid)
end


function fitdce{M,N}(Ct::Array{Float64,M}, mask::BitArray{N}, t::Vector{Float64},
                   Cp::Vector{Float64}; models=[2], residthresh=1.0, Ktmax=5.0)
  @dprint "fitting DCE data"
  sizein = size(Ct)
  n = prod(sizein[2:end])
  nt = sizein[1]
  @assert nt == length(t)
  Ct = reshape(Ct, (nt, n))
  idxs = find(mask)

  nidxs = length(idxs)
  nmodels = length(models)
  @assert nmodels > 0 "at least one model must be specified"
  resid = Inf*ones(n)
  params = zeros(3, n)
  modelmap = zeros(UInt8, n)

  if 3 in models
    @dprint "attempting Extended Tofts-Kety model"
    p0 = [0.01, 0.01, 0.01]
    f3(x,p) = extendedtoftskety(x, p, Cp)
    p, r, dof = nlsfit(f3, Ct, idxs, t, p0)
    p[2,idxs] = p[1,idxs] ./ p[2,idxs]
    resid[:] = squeeze(sum(abs2, r, 1), 1) / dof
    params[:] = p
    modelmap[idxs] = 3
  end
  if 2 in models
    @dprint "attempting Standard Tofts-Kety model"
    p0 = [0.01, 0.01]
    f2(x,p) = toftskety(x, p, Cp)
    p, r, dof = nlsfit(f2, Ct, idxs, t, p0)
    r = squeeze(sum(abs2, r, 1), 1) / dof
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
  if 1 in models
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
      params[1,k] = clamp.(params[1,k], eps(), Ktmax)
      params[2:end,k] =  clamp.(params[2:end,k], eps(), 1.0)
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



function fitdata(opts::Dict)
  @dprint "running models"

  # option order of precedence
  # 1. function arguments
  # 2. values in "datafile"
  # 3. defaults
  def = defaults()
  infile = haskey(opts,"datafile") ? opts["datafile"] : def["datafile"]
  mat = matread(infile)
  opts = merge(def, merge(mat, opts))
  validate(opts)

  # load DCE and R1 data
  relaxivity = opts["relaxivity"]
  TR = opts["TR"]
  models = opts["models"]
  outfile = opts["outfile"]
  Cp = vec(opts["Cp"])  # colon makes sure that we get a vector

  # parallel workers are better than multithreaded BLAS for this problem
  # run julia with the '-p <n>' flag to start with n workers
  BLAS.set_num_threads(1)
  # startworkers(opts["workers"])
  # require("DCEMRI.jl")

  if haskey(opts, "R10") && haskey(opts, "S0")
    @dprint "found existing R1 map"
    R10 = opts["R10"]
    S0 = opts["S0"]
  else
    @dprint "found multi-flip data"
    t1data = opts["T1data"]
    t1flip = vec(opts["T1flip"]) * pi / 180.0
    @assert length(t1flip) == size(t1data,1) # must have one flip angle per T1 image
    R10, S0 = fitr1(t1data, t1flip, TR)
  end

  dcedata = opts["DCEdata"]
  dceflip = opts["DCEflip"] * pi / 180.0
  t = vec(opts["t"]) / 60.0  # convert time to min

  dims = size(dcedata)[2:end]
  @assert dims == size(R10) "R10 map and DCE images have different dimensions"

  # MAIN: run postprocessing steps
  SER = ser(dcedata)
  if haskey(opts, "mask") # if mask specified, use it
      mask = BitArray(opts["mask"])
  else   # use SER threshold
      mask = SER .> opts["SERcutoff"]
  end
  R1 = r1eff(dcedata, R10, TR, dceflip)
  Ct = tissueconc(R1, R10, relaxivity)
  params, resid, modelmap = fitdce(Ct, mask, t, Cp, models=models)

  @dprint "saving results to $outfile"
  results = Dict()
  results["t"] = t
  results["Cp"] = Cp
  results["R10"] = R10
  results["S0"] = S0
  results["SER"] = SER
  results["mask"] = Array(mask)
  results["models"] = models
  results["modelmap"] = modelmap
  results["R1"] = R1
  results["Ct"] = Ct
  params = reshape(params, (size(params,1), prod(dims)))
  results["Kt"] = reshape(params[1,:], dims)
  results["ve"] = reshape(params[2,:], dims)
  results["vp"] = reshape(params[3,:], dims)
  results["resid"] = resid
  matwrite(outfile, results)

  results
end
fitdata(filename::AbstractString) = fitdata(datafile=filename) # point to MAT file
fitdata(; kwargs...) = fitdata(kwargs2dict(kwargs))
