function spgreqn(x::Vector{Float64}, p::Vector{Float64}, TR::Float64)
  y = similar(x)
  E = exp(-TR*p[2])  # hoisted expensive, invariant result
  @simd for k in 1:length(x)
     @inbounds y[k] = p[1]*sin(x[k])*(1.0 - E) / (1.0 - E*cos(x[k]))
  end
  y
end

function toftskety(t::Vector{Float64}, p::Vector{Float64}, Cp::Vector{Float64})
  # Standard Tofts-Kety Model.  Works when t_dce = t_aif only, so you should
  # resample the AIF to match the DCE dynamic spacing before using.
  Ct = zeros(t)
  kep = p[2]  # kep = Ktrans / ve
  for k in 1:length(Ct)
    #Ct[k] = p[1]*trapz(exp(-p[1]*(t[k] - t[1:k])/p[2]) .* Cp[1:k], t[1:k])
    # inlined trapz here ... original equations are up here ^
    tk = t[k]
    @simd for j = 1:k
      @inbounds y = exp(kep*(t[j] - tk)) * Cp[j]
      @inbounds Ct[k] += ifelse((j == 1 || j == k) && k > 1, 0.5*y, y)
    end
  end
  dtp = (t[2] - t[1]) * p[1]
  @simd for k = 1:length(Ct)
    @inbounds Ct[k] *= dtp
  end
  Ct
end

function extendedtoftskety(t::Vector{Float64}, p::Vector{Float64}, Cp::Vector{Float64})
  # Extended Tofts-Kety Model. Works when t_dce = t_aif only, so you should
  # resample the AIF to match the DCE dynamic spacing before using.
  Ct = toftskety(t, p, Cp)
  @simd for k in 1:length(Ct)
    @inbounds Ct[k] += p[3]*Cp[k]
  end
  Ct
end

function onecompartment(t::Vector{Float64}, p::Vector{Float64}, Cp::Vector{Float64})
  # plasma only model
  Ct = similar(Cp)
  @simd for k in 1:length(Ct)
    @inbounds Ct[k] = p[1]*Cp[k]
  end
  Ct
end
