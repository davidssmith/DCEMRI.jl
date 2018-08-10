function nlsfitworker(model::Function, y::Matrix{T}, x::Vector{T}, p0::Vector{T}; kwargs...) where T
  # performs Levenberg-Marquardt fitting
  nparams = length(p0)
  n = size(y,2)
  params = zeros(nparams, n)
  resids = similar(y)
  for j in 1:n
    fit = curve_fit(model, x, y[:,j], p0)
    params[:,j] = fit.param
    resids[:,j] = fit.resid
  end
  (params, resids)
end

function nlsfit(f::Function, y::Matrix{T}, idxs::Vector{Int}, x::Vector{T},
  p0::Vector{T}; kwargs...) where T
  timerStart = time_ns()
  nt, n = size(y)
  nw = nworkers()
  nidxs = length(idxs)
  params = zeros(length(p0), n)
  resids = zeros(nt, n)
  dof = nt - length(p0)
  if nw > 1
    reflist = Any[]
    idxperworker = ceil(Int, nidxs / nw)
    workeridxs = Any[idxs[(w-1)*idxperworker+1:min(w*idxperworker,nidxs)] for w in 1:nw]
    ni = map(length, workeridxs)
    @dprint "work distribution: $nt x $ni points"
    workerids = workers()
    for w in 1:nw
      r = @spawnat workerids[w] nlsfitworker(f, y[:,workeridxs[w]], x, p0; kwargs...)
      push!(reflist, r)
    end
    for w in 1:nw
      params[:,workeridxs[w]], resids[:,workeridxs[w]] = fetch(reflist[w])
    end
  else   # run on one CPU core
    @dprint "running $nt x $nidxs points on one CPU core"
    params[:,idxs], resids[:,idxs] = nlsfitworker(f, y[:,idxs], x, p0; kwargs...)
  end
  t = (time_ns()-timerStart)/1e9
  vps = nidxs / t
  @dprint @sprintf "processed %d voxels in %.1f s (%.1f vox/s)\n" nidxs t vps
  (params, resids, dof)
end

# cumtrapz() is taken from the ElasticFDA.jl package, with original license:
#
# Copyright (c) 2016: J. Derek Tucker.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 
function cumtrapz(x::Array{Float64, 1}, y::Array{Float64}, dim::Integer=1)
    # Cumulative trapezoidal integration
    # Obtained from:
    # https://github.com/jdtuck/ElasticFDA.jl/blob/master/src/misc_funcs.jl
    perm = [dim:max(length(size(y)),dim); 1:dim-1];
    y = permutedims(y, perm);
    if ndims(y) == 1
        n = 1;
        m = length(y);
    else
        m, n = size(y);
    end

    if n == 1
        dt = diff(x)/2.0;
        z = [0; cumsum(dt.*(y[1:(m-1)] + y[2:m]))];
    else
        dt = repmat(diff(x)/2.0,1,n);
        z = [zeros(1,n); cumsum(dt.*(y[1:(m-1), :] + y[2:m, :]),1)];
        z = ipermutedims(z, perm);
    end

    return z
end
