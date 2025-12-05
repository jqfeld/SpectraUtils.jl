module SpectraUtilsNonlinearSolveExt

using SpectraUtils, NonlinearSolve
using FiniteDiff, ForwardDiff, LinearAlgebra
import SpectraUtils: LevenbergMarquardtFit, fit_spectrum, FitResult, _namedtuple_statfunc
using Distributions



regularization_term(d, x) = error("Regularization term not defined for prior distribution: $(typeof(d))")
regularization_term(d::Normal, x) = (x - d.μ)/d.σ 
regularization_term(d::Uniform, x) = d.a <= x <= d.b ? 0. : Inf 


function (lm_model::LevenbergMarquardtFit)(xs, ys)
  # Priors
  if any([p isa Distribution for p in lm_model.params])
    @warn "No data uncertainty, Distributions not used as Prior."
  end

  function resid(p, _)
    pars = (; (keys(lm_model.params) .=> p)...)
    return ys .- lm_model.func(pars)(xs)
  end

  prob = NonlinearLeastSquaresProblem(
    NonlinearFunction(resid), [a for a in _namedtuple_statfunc(lm_model.params, mean)] ;)

end


function fit_spectrum(lm_model::LevenbergMarquardtFit, xs, ys;)
  prob =  lm_model(xs,ys)
  sol = solve(prob, LevenbergMarquardt(; disable_geodesic=Val(true), autodiff=AutoForwardDiff()),)
  σ² = (sum(abs2, sol.resid)/(length(xs) - length(sol.u)))
  LP(θ) = -sum(abs2, sol.prob.f(θ, [])) / 2 / σ²

  H = ForwardDiff.hessian(p -> -LP(p), sol.u)
  if any(isnan, H)
    H = FiniteDiff.finite_difference_hessian(p -> -LP(p), sol.u)
  end

  std_fit = sqrt.(abs.(diag(pinv(H))))

  return FitResult((; (keys(lm_model.params) .=> Normal.(sol.u, std_fit))...), sol.resid, nothing)
end


function (lm_model::LevenbergMarquardtFit)(xs, ys, sigmas)
  # Priors

  function resid(p, _)
    pars = (; (keys(lm_model.params) .=> p[1:end])...)
    R = [
      (ys .- lm_model.func(pars)(xs)) ./ sigmas ; 
      [regularization_term(lm_model.params[k], pars[k]) for k in keys(filter(x -> x isa Distribution, lm_model.params))]
    ]
    R
  end
  prob = NonlinearLeastSquaresProblem(
    NonlinearFunction(resid), [a for a in _namedtuple_statfunc(lm_model.params, mean)] ;)

end

function fit_spectrum(lm_model::LevenbergMarquardtFit, xs, ys, sigmas;)
  prob =  lm_model(xs,ys, sigmas)
  sol = solve(prob, LevenbergMarquardt(; disable_geodesic=Val(true), autodiff=AutoForwardDiff()),)
  
  LP(θ) = -sum(abs2, sol.prob.f(θ, [])) / 2

  H = ForwardDiff.hessian(p -> -LP(p), sol.u)
  if any(isnan, H)
    H = FiniteDiff.finite_difference_hessian(p -> -LP(p), sol.u)
  end

  std_fit = sqrt.(abs.(diag(pinv(H))))

  return FitResult((; (keys(lm_model.params) .=> Normal.(sol.u, std_fit))...), sol.resid[1:end-length(filter(x -> x isa Distribution, lm_model.params))] .* sigmas, sigmas)
end

end
