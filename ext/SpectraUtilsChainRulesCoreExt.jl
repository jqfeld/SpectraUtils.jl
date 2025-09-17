module SpectraUtilsChainRulesCoreExt

using SpectraUtils, ChainRulesCore
import SpectraUtils: gaussian, lorentzian, voigt
import SpecialFunctions: faddeeva

ChainRulesCore.@scalar_rule(gaussian(x, σ), (-x / σ^2 * gaussian(x, σ), (x^2 / σ^3 - 1 / σ) * gaussian(x, σ)))
ChainRulesCore.@scalar_rule(lorentzian(x, γ),
  (
    -2 * x * γ / π / (x^2 + γ^2)^2,
    (x^2 - γ^2) / π / (x^2 + γ^2)^2
  )
)

ChainRulesCore.@scalar_rule(voigt(x, σ, γ),
  @setup(
    z = (x + im * γ) / σ / sqrt(2),
    Im_w = imag(faddeeva(z)),
    Re_w = real(faddeeva(z)),
  ),
  (
    1 / σ^3 / sqrt(2π) * (γ * Im_w - x * Re_w),
    1 / σ^4 / sqrt(2π) * (
      (x^2 - σ^2 - γ^2) * Re_w -
      2 * x * γ * Im_w +
      γ * σ * sqrt(2 / π)
    ),
    -1 / σ^3 / sqrt(2 * π) * (σ * sqrt(2 / π) - x * Im_w - γ * Re_w),
  )
)


end
