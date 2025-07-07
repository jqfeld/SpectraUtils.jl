using SpecialFunctions

abstract type LineShape end

lineparam(x::Real, p) = x
lineparam(x::Function, p) = x(p)



struct Lorentzian{A,X,T} <: LineShape
  ampl::A
  x0::X
  hwhm::T
end

Lorentzian(x0, hwhm) = Lorentzian(1., x0, hwhm)

(L::Lorentzian{A,X,T})(x, p=[]) where {A, X, T} =  
  lineparam(L.ampl,p)/π * lineparam(L.hwhm,p) / ((x-lineparam(L.x0,p))^2 + lineparam(L.hwhm,p)^2)

struct Gaussian{A,X,T} <: LineShape
  ampl::A
  x0::X
  sigma::T
end

Gaussian(x0, sigma) = Gaussian(1., x0, sigma)

(L::Gaussian{A,X,T})(x, p=[]) where {A, X, T} = 
  lineparam(L.ampl, p)/(lineparam(L.sigma,p) * sqrt(2π)) * exp(-0.5 * ((x-lineparam(L.x0,p))/lineparam(L.sigma,p))^2)


struct VoigtApprx{A,X,T1,T2} <: LineShape
  ampl::A
  x0::X
  sigma::T1
  gamma::T2
end

VoigtApprx(x0, sigma, gamma) = VoigtApprx(1., x0, sigma, gamma)

@inline gaussian(x, σ=0.01) = 1 / sqrt(2π) / σ * exp(-x^2 / 2 / σ^2)
@inline lorenzian(x, γ) = γ / π / (γ^2 + x^2)
function (L::VoigtApprx{A,X,T})(x, p=[]) where {A, X, T} 
  μ = lineparam(L.x0, p)
  σ = lineparam(L.sigma, p)
  γ = lineparam(L.gamma, p)

  fG = 2 * sqrt(2 * log(2)) * σ
  fL = 2 * γ
  f = (max(0.0, fG^5 + 2.69269 * fG^4 * fL + 2.42843 * fG^3 * fL^2 + 4.47163 * fG^2 * fL^3 + 0.07842 * fG * fL^4 + fL^5))^(1 / 5)
  η = 1.36603 * (fL / f) - 0.47719 * (fL / f)^2 + 0.11116 * (fL / f)^3
  return lineparam(L.ampl,p)*(η * lorenzian(x - μ, f / 2) + (1 - η) * gaussian(x - μ, f / (2 * sqrt(2 * log(2)))))
end



struct Voigt{A,X,T1,T2} <: LineShape
  ampl::A
  x0::X
  sigma::T1
  gamma::T2
end

Voigt(x0, sigma, gamma) = Voigt(1., x0, sigma, gamma)

function (L::Voigt)(x, p=[])
  σ = lineparam(L.sigma, p)
  γ = lineparam(L.gamma, p)
  return real(
    faddeeva( (x + im*γ)/σ/sqrt(2) )
  ) / σ / sqrt(2π)
end

# @inline σ_doppler(f0, T; m=4.65e-26) = sqrt(k * T / m / c0^2) * f0
