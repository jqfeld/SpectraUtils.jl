module SpectraUtilsChainRulesCoreExt

using SpectraUtils, ChainRulesCore
import SpectraUtils: gaussian, lorentzian, voigt

ChainRulesCore.@scalar_rule(gaussian(x, σ), (-x/σ^2 * gaussian(x,σ), (x^2/σ^3 - 1/σ)*gaussian(x,σ)))

end
