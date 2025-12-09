using Distributions, Statistics

"""
    fit_spectrum(fit, xs, ys[, sigmas]; kwargs...)

Fit the spectral model encoded by `fit` against the observations `ys` at
positions `xs`. Dispatch is provided by extension modules: loading `Turing`
enables Bayesian inference via [`TuringFit`](@ref), while loading
`NonlinearSolve` enables deterministic fitting via
[`LevenbergMarquardtFit`](@ref). If neither extension is available, an error is
thrown indicating the missing dependency.
"""
fit_spectrum(x...) = error("Load Turing.jl to use this functionality.")


"""
    TuringFit(model, priors)

Wrap a spectral `model` and `priors` into a callable object for probabilistic
fitting with Turing. `model` should map a named tuple of parameters to a
[`Spectrum`](@ref), while `priors` holds a matching set of `Distribution`
objects. See [`fit_spectrum`](@ref) for usage examples.
"""
struct TuringFit{F,P}
  func::F
  params::P
end


"""
    LevenbergMarquardtFit(model, initial_params)

Wrap a spectral `model` and `initial_params` into a callable object for
deterministic nonlinear least-squares fitting. The `initial_params` named tuple
provides starting values (or optional priors) for each parameter used by
`model`. Used in conjunction with [`fit_spectrum`](@ref) when
`NonlinearSolve.jl` is available.
"""
struct LevenbergMarquardtFit{F,P}
  func::F
  params::P
end

"""
    FitResult(params, resid, wt)

Container returned by [`fit_spectrum`](@ref) encapsulating parameter
posteriors or estimates (`params`), residuals (`resid`), and optional weights
(`wt`). `FitResult` implements the `StatsAPI` interface, providing accessors
such as `coef`, `stderror`, `rss`, and `weights` for downstream analysis.
"""
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

