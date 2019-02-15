% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% This function counts the number of elements in vec located in bins and
% does not take into account the values outside the edges.
% Bins must be equally spaced
% As third argument, you can put 'n' to get a n normalized count (%) or 'nc' to
% get a cumulative normalized count (between 0 and 1)

function count = fct_hist(vec, bins, varargin)

if abs(std(diff(bins))) > (1000*eps)
    error('Bins are not equally spaced.');
end

pas = bins(2) - bins(1);

bins = bins + 0.5*pas;
bins = [(bins(1) - pas) bins];
bins = bins - 2*eps;

count = hist(vec, bins);
count(1) = [];
count(end) = [];


if any(strcmp(varargin, 'n'))
    count = (count / sum(count))*100;
end

if any(strcmp(varargin, 'nc'))
    count = cumsum(count);
    count = count / count(end); 
end
   
end








