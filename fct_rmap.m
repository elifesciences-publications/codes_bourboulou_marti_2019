% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% INPUT
% xpos    = animal's position over time
% idx_spk = index of spikes times
% prm     = parameters structure with fields
%    -> freq      = sample frequency
%    -> xbin      = set of positions bins
%    -> tbin      = set of times intervals
%    -> idx_rem   = index of spikes to remove (e.g., low velocity)
%    -> xbin_rem  = number of xbin to remove on both edges

% OUTPUT
% rmap = ratemap structure with matrix as fields (row: time bins (t), col: position bins (x))
%    -> nbspk_tx   = spikes count 
%    -> dwell_tx   = time spent 
%    -> fr_tx      = firing rate
%    -> nbspk_s_tx = smoothed (s) spikes count 
%    -> dwell_s_tx = smoothed (s) time spent
%    -> fr_s_tx    = smoothed (s) firing rate


function rmap = fct_rmap(xpos, idx_spk, prm)


xpos(prm.idx_rem) = NaN;
prm.xbin([1:prm.xbin_rem (end -(prm.xbin_rem - 1)):end]) = [];
idx_spk(ismember(idx_spk, prm.idx_rem)) = [];

nb_tbin = size(prm.tbin, 1);
nb_xbin = length(prm.xbin) - 1;

rmap.nbspk_tx = zeros(nb_tbin, nb_xbin);
rmap.dwell_tx = zeros(nb_tbin, nb_xbin);

for t = 1:nb_tbin
    rmap.nbspk_tx(t, :) = fct_hist(xpos(idx_spk(idx_spk >= prm.tbin(t, 1) & idx_spk <= prm.tbin(t, 2))), prm.xbin);
    rmap.dwell_tx(t, :) = (fct_hist(xpos(prm.tbin(t, 1):prm.tbin(t, 2)), prm.xbin) / prm.freq) + eps;
end

rmap.fr_tx = rmap.nbspk_tx ./ rmap.dwell_tx;


rmap.nbspk_s_tx = fct_smoothgauss(rmap.nbspk_tx, prm.ismooth);
rmap.dwell_s_tx = fct_smoothgauss(rmap.dwell_tx, prm.ismooth) + eps;
rmap.fr_s_tx = rmap.nbspk_s_tx ./ (rmap.dwell_s_tx);