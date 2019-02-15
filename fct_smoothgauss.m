% Author(s): Marti Geoffrey, Bourboulou Romain
% Epsztein Lab 2019

function data_s = fct_smoothgauss(data, hwin)

if ~ismatrix(data) && iscolumn(data)
    data = data';
end
  

data_p = padarray(data, [0 hwin], 'replicate' , 'both');


win = (-hwin:hwin);
kernel = exp(-win.^2/(hwin/2)^2);
kernel = kernel / sum(kernel);


data_s = fct_conv(data_p, kernel);
data_s = data_s(:, (2*hwin) + 1:end-(2*hwin));




