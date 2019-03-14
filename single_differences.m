function [single_diff] = single_differences(time_interval, obs_type, system, obsMatrix1, obsMatrix2, data)
load(data)
[sysPrefix, sysName] = sysGNSS(system); % dla Galileo 
    for i=1:length(obs_type)

        i_obs = find(string(obsType1(:))==obs_type(i));
        
        epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
        
        for j=1:length(epochs)
            
            i_epoch = find(cell2mat(obsTable1(:,1))==epochs(i));
            act_constellation_1 = cell2mat(obsTable1(i_epoch, i_obs+2)); 
            act_constellation_2 = cell2mat(obsTable2(i_epoch, i_obs+2));
            act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
            nsat=length(act_constellation);
            
            for k=1:nsat
                
                prn_num = act_constellation(k); % numer k-ty satelita w danej epoce
                prn_idx = prn_num - sysPrefix; % index dla kolumny w obsMatrix
                
                single_diff(j,k,i) = obsMatrix1(i_epoch, prn_idx, i_obs) - obsMatrix2(i_epoch, prn_idx, i_obs);
                
            end
            
        end
        
    end


end

