using Distributions, Statistics

fit_spectrum(x...) = error("Load Turing.jl to use this functionality.")


struct TuringFit{F,P}
  func::F
  params::P
end


struct LevenbergMarquardtFit{F,P}
  func::F
  params::P
end

struct FitResult{P,R,W}
  params::P
  resid::R
  wt::W
end

Base.getindex(fr::FitResult, k::Symbol) = fr.params[k]
_namedtuple_statfunc(nt, func) = (; [ k => func(nt[k]) for k in keys(nt)]...)

import StatsAPI

StatsAPI.coef(fr::FitResult) = _namedtuple_statfunc(fr.params, mean)
StatsAPI.coefnames(fr::FitResult) = keys(fr.params)
StatsAPI.stderror(fr::FitResult) = _namedtuple_statfunc(fr.params, Statistics.std)

StatsAPI.nobs(fr::FitResult) = length(fr.resid)
StatsAPI.dof(fr::FitResult) = nobs(fr) - length(coef(fr))
StatsAPI.rss(fr::FitResult) = sum(abs2, fr.resid)
StatsAPI.weights(fr::FitResult) = fr.wt
StatsAPI.residuals(lfr::FitResult) = lfr.resid
mse(fr::FitResult) = rss(fr) / dof(fr)
# isconverged(lsr::LsqFitResult) = lsr.converged

