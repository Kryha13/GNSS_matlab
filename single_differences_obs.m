function [single_differ] = single_differences_obs(time_interval, obs_type, system, obsMatrix1, obsMatrix2, Az_EL_R, data)
load(data)
[sysPrefix, sysName] = sysGNSS(system); % dla Galileo 
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
single_differ = [];
    for i=1:length(obs_type)

        i_obs = find(string(obsType1(:))==obs_type(i));
        
        for j=1:length(epochs)
            
            i_epoch = find(cell2mat(obsTable1(:,1))==epochs(j));
            act_constellation_1 = cell2mat(obsTable1(i_epoch, i_obs+2)); 
            act_constellation_2 = cell2mat(obsTable2(i_epoch, i_obs+2));
            act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
            nsat=length(act_constellation);
            [sat_ref, sat_idx] = find_ref_sat(Az_EL_R, j);
            sat_ref = find(act_constellation==sat_ref);
            for k=1:nsat
                
                prn_num = act_constellation(k); % numer k-ty satelita w danej epoce
                prn_idx = prn_num - sysPrefix; % index dla kolumny w obsMatrix
                
                single_diff(k,j,i) = obsMatrix1(i_epoch, prn_idx, i_obs) - obsMatrix2(i_epoch, prn_idx, i_obs);
%                 single_diff(single_diff==0) = NaN;
                % satelita referencyjny jako pierwszy
                if k==nsat
                    single_diff([1 sat_ref],j,i) = single_diff([sat_ref 1],j,i);
                end
            end
        end
        single_diff = squeeze(single_diff(:,:,i));
        single_differ = [single_differ; single_diff];   
    end
%     single_diff = reshape(single_diff.', size(single_diff, 1)* size(single_diff, 3), size(single_diff, 2)).';
end

