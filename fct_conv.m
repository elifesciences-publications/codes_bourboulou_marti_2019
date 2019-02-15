% Author(s): Marti Geoffrey, Bourboulou Romain
% Epsztein Lab 2019

function out = fct_conv(x, y)

N = size(x, 2) + size(y, 2) - 1;
N2 = pow2(nextpow2(N));
 
x_m = isnan(x);
y_m = isnan(y);
x(x_m) = 0;
y(y_m) = 0;

x_f = fft(x, N2 , 2); 
y_f = fft(y, N2 , 2);

out = ifft(bsxfun(@times , x_f, y_f), N2, 2);
out = out(:, 1:N);
out(x_m) = NaN;
