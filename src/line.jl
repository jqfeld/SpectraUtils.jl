"""
    Line(amplitude, position, shape)

Bundle a spectral line definition with an `amplitude`, `position`, and callable
`shape`. Each field may be either a numeric value or a function that consumes a
parameter set when the line is evaluated. Call the resulting object as
`line(x, p)` to evaluate the profile at a scalar position and as `line(xs, p)` to
broadcast over a vector of positions.
"""
struct Line{A,P,LS}
  amplitude::A
  position::P
  shape::LS
end

@inline (l::Line{A,P,LS})(x, p=NullParameters()) where {A,P,LS} =
  calc_param(l.amplitude, p) * l.shape(x - calc_param(l.position, p), p)

function (l::Line{A,P,LS})(xs::Vector{X}, p=NullParameters()) where {A,P,LS,X}
  out = zeros(typeof(l(xs[1], p)), length(xs))
  @simd for i in eachindex(xs)
    out[i] = calc_param(l.amplitude, p) * l.shape(xs[i] - calc_param(l.position, p), p)
  end
  return out
end
