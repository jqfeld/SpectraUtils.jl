module SpectraUtilsTuringExt

using SpectraUtils, Turing
import SpectraUtils: TuringFit, fit_spectrum, FitResult
using Distributions


@model function (turing_model::TuringFit)(xs, ys)
  # Priors
  params ~ product_distribution(turing_model.params)
  sigma ~ turing_model.sigma

  # Likelihood
  spec = turing_model.func(params)
  μ = spec.(xs)

  ys ~ MvNormal(μ, sigma^2 * I)
end



function fit_spectrum(turing_model::TuringFit, xs, ys;
  sampler=NUTS(), num_samples=1000, num_chains=1, progress=true)

  posterior = turing_model(xs, ys)
  chain = sample(posterior, sampler, MCMCThreads(), num_samples, num_chains; progress)
  ks = Symbol.(replace.(namesingroup(chain, "params", index_type=:dot) .|> string, "params." => ""))
  return FitResult((; [k => vec(chain[Symbol("params.", k)]) for k in ks]...), vec(chain[:sigma]), nothing)
end

end
