using Distributions, Statistics
# import StatsAPI

fit_spectrum(x...) = error("Load Turing.jl to use this functionality.")


struct TuringFit{F,P,S}
  func::F
  params::P
  sigma::S
end

TuringFit(func, params) = TuringFit(func,params, Gamma(1,0.5))

struct FitResult{P,S,R}
  params::P
  sigma::S
  resid::R
end

_namedtuple_statfunc(nt, func) = (; [ k => func(nt[k]) for k in keys(nt)]...)

coef(fr::FitResult) = _namedtuple_statfunc(fr.params, mean)
coefnames(fr::FitResult) = keys(fr.params)
stderror(fr::FitResult) = _namedtuple_statfunc(fr.params, Statistics.std)

Base.getindex(fr::FitResult, k::Symbol) = fr.params[k]


# StatsAPI.dof(fr::FitResult) = nobs(fr) - length(coef(lfr))
# StatsAPI.nobs(fr::FitResult) = length(fr.resid)
# StatsAPI.rss(fr::FitResult) = sum(abs2, fr.resid)
# StatsAPI.weights(lfr::FitResult) = fr.wt
# StatsAPI.residuals(lfr::FitResult) = lfr.resid
# mse(lfr::LsqFitResult) = rss(lfr) / dof(lfr)
# isconverged(lsr::LsqFitResult) = lsr.converged

