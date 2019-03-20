function [single_differ] = single_diff_u(time_interval, obs1, obs2, Az_EL_R, d_wsp, data)
load(data)
% [sysPrefix, sysName] = sysGNSS(system); % dla Galileo 
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
single_diff = {};
single_differ = [];
j = 1;
    for i=1:length(epochs)
    
        act_constellation_1 = squeeze(ro_approx1(:,2,i)); 
        act_constellation_2 = squeeze(ro_approx2(:,2,i));
        act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
        nsat=length(act_constellation);
        [sat_ref, sat_idx] = find_ref_sat(Az_EL_R, i); % satelita referencyjny
        % tu moze byc problem jak cos 
        sat_ref = find(act_constellation==sat_ref); % indeks w aktualnej konstelacji
        
        for k=1:nsat
                
            prn_num = act_constellation(k); % numer k-ty satelita w danej epoce
            idx_1 = find(squeeze(ro_approx1(:,2,i))==prn_num);
            idx_2 = find(squeeze(ro_approx2(:,2,i))==prn_num);

            single_diff(k,i) = {obs1(idx_1, d_wsp, i) - obs2(idx_2, d_wsp, i)};
            
            if (single_diff{k,i}==0)
                single_differ(k,i) = 0;
            else
                single_differ(k,i) = squeeze(single_diff{k,i});
            end
            % satelita referencyjny jako pierwszy
            if k==nsat
                single_differ([1 sat_ref],i) = single_differ([sat_ref 1],i);
            end
        end       
    end
end

