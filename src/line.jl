"""
    Line(amplitude, position, shape)

Bundle a spectral line definition with an `amplitude`, `position`, and callable
`shape`. Each field may be either a numeric value or a function that consumes a
parameter set when the line is evaluated. Call the resulting object as
`line(x, p)` to evaluate the profile at a scalar position and as `line(xs, p)` to
broadcast over a vector of positions.
"""
struct Line{A,P,LS}
  position::P
  amplitude::A
  shape::LS
end


@inline (l::Line{A,P,LS})(x) where {A,P,LS} = l.amplitude * l.shape(x - l.position)

function (l::Line{A,P,LS})(xs::Vector{X}) where {A,P,LS,X}
  out = zeros(typeof(l(xs[1])), length(xs))
  @simd for i in eachindex(xs)
    out[i] = l.amplitude * l.shape(xs[i] - l.position)
  end
  return out
end
