% Author(s): Marti Geoffrey
% Epsztein Lab 2019


function ivec = fct_discretize(vec, bins, varargin)


[vec, bins, prm] = dealin(vec, bins, varargin);


% vec = [-10 -15 -14 0.1 0.2 0.3 2.1 2.2 2.3 2.56 3.1 4 10 15 20 30 12];
% bins = [0 1 2 3 4];

% vec = [-5 0 0.1 0.8 1.5 1.8 3.1 3.5 4.5 10 15 16 18 25 59 70 78];
% bins = [0 1 ; 3 9 ; 20 60];


switch prm.mtd
    case 'closest' % look for the closest values from the vector vec in the vector bins
        switch prm.algo
            case 1 % old algo, but with "isedgenan" option, values have to be equally spaced
                if prm.isedgenan % NaN for values outside the edges
                    ts = mean(diff(bins));
                    bins = [(bins(1) - ts) bins];
                    bins = [bins (bins(end) + ts)];
                    
                    ind = 0:1:(length(bins) - 1);
                    ind(1) = NaN;
                    ind(end) = NaN;
                else
                    ind = 1:length(bins);
                end
                
                [vec_perm, perm] = sort(vec);
                kao = hist(vec_perm, bins);
                ind2 = find(kao >= 1);
                ivec = [];
                
                if ~isempty(ind2)
                    for k = ind2
                        ivec = [ivec repmat(ind(k), 1, kao(k))];
                    end
                end
                
                ivec = ivec(fct_inv_perm(perm));
                
            case 2 % values have to be equally spaced
                ts = mean(diff(bins));
                bins = bins - 0.5*ts;
                bins = [bins (bins(end) + ts)];
                
                bins(1) = -inf;
                bins(end) = +inf;
                
                [~, ivec] = histc(vec, bins);
                
                
            case 3 % the most appropriate since it's intented to find the closest values by error minimization but needs much memory
                
                tmp = abs(bsxfun(@minus, vec, bins'));
                [~, ivec] = min(tmp);   
        end
        
               
    case 'bin'
        
        if prm.isedgenan
        else
            bins(1) = -inf;
            bins(end) = +inf;
        end
        
        [~, ivec] = histc(vec, bins);
        ivec(ivec == 0) = NaN;
        
    case 'disjint'
        
        if ~prm.nointernan
            [~, ivec] = histc(vec, bins);
            
            ivec(mod(ivec, 2) == 0) = NaN; % NaN for unlabeled intervals
            ivec(mod(ivec, 2) == 1) = floor(ivec(mod(ivec, 2) == 1) / 2) + 1;
            
        else
            bins(1) = -inf;
            bins(end) = +inf;
            
            [~, ivec] = histc(vec, bins);
            
            ivec(mod(ivec, 2) == 0) = ivec(mod(ivec, 2) == 0) - 1; % Approximation...
            ivec(mod(ivec, 2) == 1) = floor(ivec(mod(ivec, 2) == 1) / 2) + 1;
        end
        
        
    case 'interint'
        warning('When the input is a set of intervals, they have to be disjoints.')
        ivec = NaN(1, length(vec));
end

ivec = ivec';

end


function [vec, bins, prm] = dealin(vec, bins, X)

if size(vec, 1) > 1
    vec = vec';
end


if any(strcmp(X, 'nointernan'))
    prm.nointernan = true;
else
    prm.nointernan = false;
end


if size(bins, 1) > 1 && size(bins, 2) > 1
    bins = bins';
    bins = bins(:)';
    diffbins = diff(bins);
    if all(diffbins(2:2:end) > 0) % Disjoints Intervals
        prm.mtd = 'disjint';
        return
    else
        prm.mtd = 'interint';
    end
end


if size(bins, 1) > 1
    bins = bins';
end


if any(strcmp(X, 'edgenan'))
    prm.isedgenan = true;
    prm.algo = 1;
else
    prm.isedgenan = false;
    if (abs(std(diff(bins))) <= 1000*eps)
        prm.algo = 2;
    else
        warning('Bins are not equally spaced')
        prm.algo = 3;
    end
end


if any(strcmp(X, 'bin'))
    prm.mtd = 'bin';
else
    prm.mtd = 'closest';
end






end


