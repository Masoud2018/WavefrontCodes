function g = Fnear_comp(f0,H)

g = ifourier2D(fourier2D(f0) .* H);

end