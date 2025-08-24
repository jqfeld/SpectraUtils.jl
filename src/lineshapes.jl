using SpecialFunctions

abstract type LineShape end

@inline lorentzian(x, γ) = γ / π / (γ^2 + x^2)

struct Lorentzian{T} <: LineShape
  hwhm::T
end


(L::Lorentzian{T})(x, p=[]) where {T} = lorentzian(x, calc_param(L.hwhm, p))  

@inline gaussian(x, σ) = 1 / sqrt(2π) / σ * exp(-x^2 / 2 / σ^2)

struct Gaussian{T} <: LineShape
  sigma::T
end

(L::Gaussian{T})(x, p=[]) where {T} = gaussian(x, calc_param(L.sigma,p))


struct VoigtApprx{A,X,T1,T2} <: LineShape
  sigma::T1
  gamma::T2
end

VoigtApprx(sigma::T1, gamma::T2) where {T1,T2} = VoigtApprx{nothing,nothing,T1,T2}(sigma, gamma)


function (L::VoigtApprx)(x, p=[])
  σ = calc_param(L.sigma, p)
  γ = calc_param(L.gamma, p)

  fG = 2 * sqrt(2 * log(2)) * σ
  fL = 2 * γ
  f = (max(0.0, fG^5 + 2.69269 * fG^4 * fL + 2.42843 * fG^3 * fL^2 + 4.47163 * fG^2 * fL^3 + 0.07842 * fG * fL^4 + fL^5))^(1 / 5)
  η = 1.36603 * (fL / f) - 0.47719 * (fL / f)^2 + 0.11116 * (fL / f)^3
  return η * lorentzian(x, f / 2) + (1 - η) * gaussian(x, f / (2 * sqrt(2 * log(2))))
end



struct Voigt{T1,T2} <: LineShape
  sigma::T1
  gamma::T2
end


voigt(x, σ, γ) = real(
    faddeeva( (x + im*γ)/σ/sqrt(2) )
  ) / σ / sqrt(2π)

function (L::Voigt)(x, p=[])
  σ = calc_param(L.sigma, p)
  γ = calc_param(L.gamma, p)
  return voigt(x, σ, γ)
end

