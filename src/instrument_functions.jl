using QuadGK

gaussian(x, sigma) = 1 / sigma / sqrt(2π) * exp(-0.5 * (x / sigma)^2)
struct GaussianInstrument
  sigma
end
(g::GaussianInstrument)(x) = gaussian(x, g.sigma)
Base.length(::GaussianInstrument) = 1
Base.iterate(g::GaussianInstrument) = (g, nothing) # Dummy iterator for compatibility
Base.iterate(g::GaussianInstrument, i) = nothing # Dummy iterator for compatibility


function convolve_instrument(x, inst::GaussianInstrument, line::Voigt)
  linewidth = max(line.sigma, line.gamma)
  if abs(x - line.x0) > max(inst.sigma, max(line.sigma, line.gamma)) * 10 
    return 0.
  else
    return quadgk(y -> inst(y+(line.x0-x))*line(y+line.x0), -linewidth*10, linewidth*10, rtol=1e-3)[1]
  end
end
