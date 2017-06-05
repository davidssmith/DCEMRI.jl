
function levenberg_marquardt{T}( f::Function, g::Function, initial_x::AbstractVector{T};
    tolX::Real = 1e-8, tolG::Real = 1e-12, maxIter::Integer = 100,
    lambda::Real = 10.0, lambda_increase::Real = 10., lambda_decrease::Real = 0.1,
    min_step_quality::Real = 1e-3, good_step_quality::Real = 0.75,
    show_trace::Bool = false, lower::Vector{T} = Array{T}(0), upper::Vector{T} = Array{T}(0))
	# finds argmin sum(f(x).^2) using the Levenberg-Marquardt algorithm
	#          x
	# The function f should take an input vector of length n and return an output vector of length m
	# The function g is the Jacobian of f, and should be an m x n matrix
	# initial_x is an initial guess for the solution
	# fargs is a tuple of additional arguments to pass to f
	# available options:
	#   tolX - search tolerance in x
	#   tolG - search tolerance in gradient
	#   maxIter - maximum number of iterations
	#   lambda - (inverse of) initial trust region radius
	#   show_trace - print a status summary on each iteration if true
	# returns: x, J
	#   x - least squares solution for x
	#   J - estimate of the Jacobian of f at x
    # check parameters
    ((isempty(lower) || length(lower)==length(initial_x)) && (isempty(upper) || length(upper)==length(initial_x))) ||
            throw(ArgumentError("Bounds must either be empty or of the same length as the number of parameters."))
    ((isempty(lower) || all(initial_x .>= lower)) && (isempty(upper) || all(initial_x .<= upper))) ||
            throw(ArgumentError("Initial guess must be within bounds."))
    (0 <= min_step_quality < 1) || throw(ArgumentError(" 0 <= min_step_quality < 1 must hold."))
    (0 < good_step_quality <= 1) || throw(ArgumentError(" 0 < good_step_quality <= 1 must hold."))
    (min_step_quality < good_step_quality) || throw(ArgumentError("min_step_quality < good_step_quality must hold."))


	const MAX_LAMBDA = 1e16 # minimum trust region radius
	const MIN_LAMBDA = 1e-16 # maximum trust region radius
	const MIN_DIAGONAL = 1e-6 # lower bound on values of diagonal matrix used to regularize the trust region step
	converged = false
	iterCt = 0
	x = initial_x
	delta_x = copy(initial_x)
	f_calls = 0
	g_calls = 0
	fcur = f(x)
	f_calls += 1
	residual = norm(fcur)^2

    # Create buffers
    n = length(x)
    JJ = Matrix{T}(n, n)
    n_buffer = Vector{T}(n)

	while ~converged && iterCt < maxIter
		J = g(x)
		g_calls += 1
		# we want to solve:
		#    argmin 0.5*||J(x)*delta_x + f(x)||^2 + lambda*||diagm(J'*J)*delta_x||^2
		# Solving for the minimum gives:
		#    (J'*J + lambda*DtD) * delta_x == -J^T * f(x), where DtD = diagm(sum(J.^2,1))
		# Where we have used the equivalence: diagm(J'*J) = diagm(sum(J.^2, 1))
		# It is additionally useful to bound the elements of DtD below to help
		# prevent "parameter evaporation".

		# REPLACED: DtD = diagm(vec(Float64[max(x, MIN_DIAGONAL) for x in sum(J.^2,1)]))
        DtD = vec(sum(abs2, J, 1))
        for i in 1:length(DtD)
          if DtD[i] <= MIN_DIAGONAL
            DtD[i] = MIN_DIAGONAL
          end
        end

		# REPLACED: delta_x = ( J'*J + sqrt(lambda)*DtD ) \ -J'*fcur
        At_mul_B!(JJ, J, J)
        @simd for i in 1:n
          @inbounds JJ[i, i] += lambda * DtD[i]
        end
        At_mul_B!(n_buffer, J, fcur)
        scale!(n_buffer, -1)
        delta_x = JJ \ n_buffer

        # apply box constraints
        if !isempty(lower)
          @simd for i in 1:n
            @inbounds delta_x[i] = max(x[i] + delta_x[i], lower[i]) - x[i]
          end
        end
        if !isempty(upper)
            @simd for i in 1:n
               @inbounds delta_x[i] = min(x[i] + delta_x[i], upper[i]) - x[i]
            end
        end

		# if the linear assumption is valid, our new residual should be:
		# REPLACED: predicted_residual = norm(J*delta_x + fcur)^2
        predicted_residual = sum(abs2, J*delta_x + fcur)

		# check for numerical problems in solving for delta_x by ensuring that the predicted residual is smaller
		# than the current residual
		if predicted_residual > residual + 2max(eps(predicted_residual),eps(residual))
			warn("""Problem solving for delta_x: predicted residual increase.
                             $predicted_residual (predicted_residual) >
                             $residual (residual) + $(eps(predicted_residual)) (eps)""")
		end
		# try the step and compute its quality
		trial_f = f(x + delta_x)
		f_calls += 1
		trial_residual = norm(trial_f)^2
		# step quality = residual change / predicted residual change
		rho = (trial_residual - residual) / (predicted_residual - residual)
		if rho > min_step_quality
			x += delta_x
			fcur = trial_f
			residual = trial_residual
			# increase trust region radius
			lambda = max(0.1*lambda, MIN_LAMBDA)
		else
			# decrease trust region radius
			lambda = min(10*lambda, MAX_LAMBDA)
		end
		iterCt += 1
		# check convergence criteria:
		# 1. Small gradient: norm(J^T * fcur, Inf) < tolG
		# 2. Small step size: norm(delta_x) < tolX
		converged = (norm(J' * fcur, Inf) < tolG) || (norm(delta_x) < tolX*(tolX + norm(x)))
	end
	x
end

function nlsfitworker{T}(model::Function, y::Matrix{T}, x::Vector{T},
  p0::Vector{T}, plower::Vector{T}=Array{T}(0), pupper::Vector{T}=Array{T}(0))
  # performs Levenberg-Marquardt fitting
  nparams = length(p0)
  n = size(y,2)
  params = zeros(nparams, n)
  resids = similar(y)
  for j in 1:n
    f(p) = model(x, p) - y[:,j]
    g = Calculus.jacobian(f, :forward)
    results = levenberg_marquardt(f, g, p0; tolX=1e-3, tolG=1e-4,
                 lambda=10.0, upper=pupper, lower=plower)
    resids[:,j] = f(results) / sqrt(norm(y[:,j]))
    params[:,j] = results
  end
  (params, resids)
end

function nlsfit{T}(f::Function, y::Matrix{T}, idxs::Vector{Int}, x::Vector{T},
  p0::Vector{T}; plower::Vector{T}=Array{T}(0), pupper::Vector{T}=Array{T}(0))
  tic()
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
      r = @spawnat workerids[w] nlsfitworker(f, y[:,workeridxs[w]], x, p0, plower, pupper)
      push!(reflist, r)
    end
    for w in 1:nw
      params[:,workeridxs[w]], resids[:,workeridxs[w]] = fetch(reflist[w])
    end
  else   # run on one CPU core
    @dprint "running $nt x $nidxs points on one CPU core"
    params[:,idxs], resids[:,idxs] = nlsfitworker(f, y[:,idxs], x, p0, plower, pupper)
  end
  t = toq()
  vps = nidxs / t
  @dprint @sprintf "processed %d voxels in %.1f s (%.1f vox/s)\n" nidxs t vps
  (params, resids, dof)
end
