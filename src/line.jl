
struct Line{A,P,LS}
  amplitude::A
  position::P
  shape::LS
end



@inline (l::Line{A,P,LS})(x, p=NullParameters()) where {A,P,LS} = calc_param(l.amplitude, p) * l.shape(x - calc_param(l.position, p), p)

function (l::Line{A,P,LS})(xs::Vector{X}, p=NullParameters()) where {A,P,LS,X}
  out = zeros(typeof(l(xs[1], p)), length(xs))
  @simd for i in eachindex(xs)
    out[i] = calc_param(l.amplitude, p) * l.shape(xs[i] - calc_param(l.position, p), p)
  end
  return out
end

