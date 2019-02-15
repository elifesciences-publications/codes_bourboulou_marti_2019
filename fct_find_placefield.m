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
%    -> nb_rep    = number of simulated spike trains

% OUTPUT
% rmap = ratemap structure with matrix as fields (row: time bins (t), col: position bins (x))
%    -> nbspk_tx   = spikes count 
%    -> dwell_tx   = time spent 
%    -> fr_tx      = firing rate
%    -> nbspk_s_tx = smoothed (s) spikes count 
%    -> dwell_s_tx = smoothed (s) time spent
%    -> fr_s_tx    = smoothed (s) firing rate
%    -> pval_pf_tx = p-val for placefield detection

function rmap = fct_find_placefield(xpos, idx_spk, prm)

rmap = fct_rmap(xpos, idx_spk, prm);
rmap.pval_pf_tx = fct_placefield_pval(xpos, rmap, prm);
