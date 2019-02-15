% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% INPUT
% xpos    = animal's position over time
% rmap    = rate map structure from "fct_rmap" function
% prm     = must be the same parameters used for the "fct_rmap" function
%    -> nb_rep    = number of simulated spike trains

% OUTPUT
%   pval_pf_tx = p-val for placefield detection (row: time bins (t), col: position bins (x))
%   fr_s_txrep = 3D matrix corresponding to a set of simulated rate maps (row: time bins (t), col: position bins (x), depth: ratemap number)


function [pval_pf_tx, fr_s_txrep] = fct_placefield_pval(xpos, rmap, prm)


xpos(prm.idx_rem) = NaN;
prm.xbin([1:prm.xbin_rem (end -(prm.xbin_rem - 1)):end]) = [];

nb_tbin = size(prm.tbin, 1);
nb_xbin = length(prm.xbin) - 1;
nb_rep = prm.nb_rep;


dwell_s_tx = rmap.dwell_s_tx;
nbspk_tx = rmap.nbspk_tx;
fr_s_tx = rmap.fr_s_tx;
tbin = prm.tbin;
xbin = prm.xbin;
ismooth = prm.ismooth;

fr_s_txrep = zeros(nb_tbin, nb_xbin, nb_rep);
pval_pf_tx = ones(nb_tbin, nb_xbin);




parfor t = 1:nb_tbin   
    idx = tbin(t, 1):tbin(t, 2);
    idxnonan = idx(1) + find(~isnan(xpos(idx))) - 1;
    N = length(idxnonan);
    nbspk = sum(nbspk_tx(t, :));
    
    nbspk_trep = poissrnd(nbspk / N, [nb_rep N]);
 
    idx_xbin = fct_discretize(xpos(idxnonan), xbin, 'bin');
    idx_xbin = repmat(idx_xbin', nb_rep, 1);
    idx_rep = repmat((1:nb_rep)', 1, N);
    
    nbspk_xrep = accumarray([idx_rep(:) idx_xbin(:)], nbspk_trep(:), [nb_rep nb_xbin]);  
    nbspk_s_xrep = fct_smoothgauss(nbspk_xrep, ismooth);
    
    fr_xrep = bsxfun(@rdivide, nbspk_s_xrep, dwell_s_tx(t, :));
    fr_s_xrep = fct_smoothgauss(fr_xrep, ismooth)';
    
    fr_s_txrep(t, :, :) = fr_s_xrep;
    pval_pf_tx(t, :) =  sum(bsxfun(@ge, fr_s_xrep, fr_s_tx(t, :)'), 2) / nb_rep;
end