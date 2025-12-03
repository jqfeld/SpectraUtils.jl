module SpectraUtilsTuringExt

using SpectraUtils, Turing
import SpectraUtils: TuringFit, fit_spectrum, FitResult, _namedtuple_statfunc

using Distributions

# Fitting spectra with unknown uncertainties

@model function (turing_model::TuringFit)(xs, ys, sigma::S=Gamma(1, 0.5)) where {S<:Distribution}
  # Priors
  params ~ product_distribution(turing_model.params)
  σ ~ sigma

  # Likelihood
  spec = turing_model.func(params)
  μ = spec.(xs)

  ys ~ MvNormal(μ, σ^2 * I)
end


function fit_spectrum(turing_model::TuringFit, xs, ys, sigma::S=Gamma(1, 0.5);
  sampler=NUTS(), num_samples=1000, num_chains=1, progress=true) where {S<:Distribution}

  posterior = turing_model(xs, ys, sigma)
  chain = sample(posterior, sampler, MCMCThreads(), num_samples, num_chains; progress)
  ks = Symbol.(replace.(namesingroup(chain, "params", index_type=:dot) .|> string, "params." => ""))
  ps = (; [k => vec(chain[Symbol("params.", k)]) for k in ks]...)
  return FitResult(ps,
    turing_model.func(_namedtuple_statfunc(ps, mean))(xs) .- ys,
    nothing
  )
end


# Fitting spectra with experimental uncertainties


@model function (turing_model::TuringFit)(xs, ys, sigmas::S) where {S<:AbstractArray}
  # Priors
  params ~ product_distribution(turing_model.params)

  # Likelihood
  spec = turing_model.func(params)
  μ = spec.(xs)

  ys ~ MvNormal(μ, diagm(sigmas .^ 2))
end


function fit_spectrum(turing_model::TuringFit, xs, ys, sigmas::S;
  sampler=NUTS(), num_samples=1000, num_chains=1, progress=true) where {S<:AbstractArray}

  posterior = turing_model(xs, ys, sigmas)
  chain = sample(posterior, sampler, MCMCThreads(), num_samples, num_chains; progress)
  ks = Symbol.(replace.(namesingroup(chain, "params", index_type=:dot) .|> string, "params." => ""))
  ps = (; [k => vec(chain[Symbol("params.", k)]) for k in ks]...)
  return FitResult(ps,
    turing_model.func(_namedtuple_statfunc(ps, mean))(xs) - ys,
    sigmas
  )
end



end
