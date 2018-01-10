function f = lanczos(x,a)
% See Graphics Gems, Andrew S. Glasser (ed), Morgan Kaufman, 1990,
% pp. 156-158.

f = (sin(pi*x) .* sin(pi*x/a) + eps) ./ ((pi^2 * x.^2 / a) + eps);
f = f .* (abs(x) < a);

end