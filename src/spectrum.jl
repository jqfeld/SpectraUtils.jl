struct Spectrum{T,B}
  lines::T
  background::B
end

Spectrum(lines) = Spectrum(lines, nothing)

function (spec::Spectrum)(x::X) where X<:Real
  ret = spec.background(x)
  for line in spec.lines
    ret += line(x)
  end
  return ret
end

function (spec::Spectrum)(x)
  ret = spec.background.(x)
  for line in spec.lines
    ret .+= line.(x)
  end
  return ret
end


function (spec::Spectrum{T,Nothing})(x::X) where {T, X<:Real}
  ret = spec.lines[1](x)
  for line in spec.lines[2:end]
    ret += line(x)
  end
  return ret
end

function (spec::Spectrum{T,Nothing})(x) where T
  ret = spec.lines[1].(x)
  for line in spec.lines[2:end]
    ret .+= line.(x)
  end
  return ret
end
