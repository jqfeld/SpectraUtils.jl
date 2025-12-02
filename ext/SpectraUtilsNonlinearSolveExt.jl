module SpectraUtilsNonlinearSolveExt

using SpectraUtils, NonlinearSolve
using ForwardDiff, LinearAlgebra
import SpectraUtils: LevenbergMarquardtFit, fit_spectrum, FitResult, _namedtuple_statfunc
using Distributions



regularization_term(d, x) = error("Regularization term not defined for prior distribution: $(typeof(d))")
regularization_term(d::Normal, x) = (x - d.μ)/d.σ 
regularization_term(d::Uniform, x) = d.a <= x <= d.b ? 0. : Inf 

function (lm_model::LevenbergMarquardtFit)(xs, ys)
  # Priors

  function resid(p, _)
    pars = (; (keys(lm_model.params) .=> p[1:end-1])...)
    σ = p[end]
    R = [
      (ys .- lm_model.func(pars)(xs)) ./ σ; 
      [regularization_term(lm_model.params[k], pars[k]) for k in keys(filter(x -> x isa Distribution, lm_model.params))];
      regularization_term(lm_model.sigma, σ)
    ]
    R
  end
  prob = NonlinearLeastSquaresProblem(
    NonlinearFunction(resid), [[a for a in _namedtuple_statfunc(lm_model.params, mean)]; mean(lm_model.sigma)] ;)

end



function fit_spectrum(lm_model::LevenbergMarquardtFit, xs, ys;)
  prob =  lm_model(xs,ys)
  sol = solve(prob, LevenbergMarquardt(; disable_geodesic=Val(true), autodiff=AutoFiniteDiff()),)
  
  LP(θ) = -sum(abs2, sol.prob.f(θ, [])) / 2

  H = ForwardDiff.hessian(p -> -LP(p), sol.u)
  if any(isnan, H)
    H = FiniteDiff.finite_difference_hessian(p -> -LP(p), sol.u)
  end

  std_fit = sqrt.(abs.(diag(pinv(H))))

  return FitResult((; (keys(lm_model.params) .=> Normal.(sol.u[1:end-1], std_fit[1:end-1]))...), Normal(sol.u[end], std_fit[end]), nothing)
end

end
