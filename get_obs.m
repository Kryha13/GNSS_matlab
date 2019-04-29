function [CLa, CLb] = get_obs(time_interval, obs_types) 

load('dane_do_test.mat')
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
fail_sats = [20 22]; % ew 14 i 18 

CLa = nan(length(epochs), 40, 8);
CLb = nan(length(epochs), 40, 8);

    for i=1:length(epochs)
        
        i_epoch = find(cell2mat(obsTable1(:,1))==epochs(i)); % indeks wiersza w obsMatrix
        
        for x = 1:14
       
            con1 = find(~isnan(obsMatrix1(i_epoch, :, x)));
            cons1 = find(~isnan(obsMatrix1(i_epoch, :, x+1)));
            constellation1{x} = intersect(con1, cons1, 'stable');
            con2 = find(~isnan(obsMatrix2(i_epoch, :, x)));
            cons2 = find(~isnan(obsMatrix2(i_epoch, :, x+1)));
            constellation2{x} = intersect(con2, cons2, 'stable');

        end
    
        for cc = 1:14

            constel{cc} = intersect(cell2mat(constellation1(cc)), cell2mat(constellation2(cc)), 'stable');

        end

    
        act_constellation = mintersect(cell2mat(constel(1)), cell2mat(constel(2)), cell2mat(constel(3)),...
            cell2mat(constel(4)),cell2mat(constel(5)),cell2mat(constel(6)),cell2mat(constel(7)),...
            cell2mat(constel(8)), cell2mat(constel(9)), cell2mat(constel(10)), cell2mat(constel(11)),...
            cell2mat(constel(12)), cell2mat(constel(13)), cell2mat(constel(14)));
        act_constellation = setdiff(act_constellation, fail_sats); % usuniêcie fai
        
        for j=1:length(obs_types)

            i_obs = find(string(obsType1(:))==obs_types(j)); % indeks warstwy w obsMatrix
            
            obs1 = obsMatrix1(i_epoch, act_constellation, i_obs);
            obs2 = obsMatrix2(i_epoch, act_constellation, i_obs);
            
            CLa(i, act_constellation, j) = obs1;
            CLb(i, act_constellation, j) = obs2; 
        end          
    end
end


