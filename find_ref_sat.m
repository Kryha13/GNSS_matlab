function [sat_ref, sat_idx] = find_ref_sat(Az_EL_R, epoch_idx)
%%%%%% INPUT %%%%%%
% Az_EL_R - tablica z azymutem, elewacja i pseudoodleg³oœci¹ dla kolejnych
% epok
% epoch_idx - interesuj¹ca nas epoka, dla której mamy wybraæ satelitê
% referencyjnego
%%%%%% OUTPUT %%%%%%
% sat_ref - numer satelity referencyjnego
% sat_idx - index satelity referencyjnego w konstelacji epoki - wiersz
%%
el = Az_EL_R(:,4,epoch_idx);
max_el = max(el);
sat_idx = find(el==max_el);
sat_ref = Az_EL_R(sat_idx,2,epoch_idx);

end

