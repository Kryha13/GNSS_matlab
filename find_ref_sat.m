function [sat_idx] = find_ref_sat(sats_elev)
%%%%%% INPUT %%%%%%
% sats_elev - wysokoœci satelit obserwowanych w danej epoce
%%%%%% OUTPUT %%%%%%
% sat_idx - index satelity referencyjnego w konstelacji epoki
%%
max_el = max(sats_elev);
sat_idx = find(sats_elev==max_el);

end

