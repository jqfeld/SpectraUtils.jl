
struct Line{A,P,LS}
  amplitude::A
  position::P
  shape::LS
end



(l::Line)(x, p) = calc_param(l.amplitude, p) *l.shape(x - calc_param(l.position,p), p)
