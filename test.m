const;
% time_interval = 100;
X1 = [3655333.1790; 1403901.2258; 5018038.3565]; %SEPT
X2 = [3655336.8160; 1403899.0382; 5018036.4527]; % SEP2
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
act_time1 = epochs + dtrec1(:,2);
act_time2 = epochs + dtrec2(:,2);
fodw = 298.2572221;
tau = 0.07; % przybli¿ony 
obs_types = ["C1C" "L1C" "C6C" "L6C" "C5Q" "L5Q" "C7Q" "L7Q" "C8Q" "L8Q"]; % wybór obserwacji 
phase_freq = [fE1 fE6 fE5 fE7 fE8];
system = 'E';
[sysPrefix, sysName] = sysGNSS(system);
fail_sats = [220 222]; % odbierany by³ sygna³ z 220 w obs a w efemerydach go nie ma nawet

[B1, L1, H1] = togeod(a, fodw, X1(1), X1(2), X1(3));
[B2, L2, H2] = togeod(a, fodw, X2(1), X2(2), X2(3));

% wybraæ tylko te które s¹ w obs_types
C1a = []; 
C1b = [];
L1a = [];
L1b = [];
C6a = [];
C6b = [];
L6a = [];
L6b = [];
C5a = [];
C5b = [];
L5a = [];
L5b = [];
C7a = [];
C7b = [];
L7a = [];
L7b = [];
C8a = [];
C8b = [];
L8a = [];
L8b = [];

del_1 = [];
del_2 = [];

for i=1:length(epochs)
    
    i_epoch = find(cell2mat(obsTable1(:,1))==epochs(i)); % indeks wiersza w obsMatrix
    
    for j=1:length(obs_types)
        
        i_obs = find(string(obsType1(:))==obs_types(j)); % indeks warstwy w obsMatrix
        % konstelacje dla obu odbiorników w danej epoce
        act_constellation_1 = cell2mat(obsTable1(i_epoch, i_obs+2)); 
        act_constellation_2 = cell2mat(obsTable2(i_epoch, i_obs+2));
        % czêœæ wspólna konstelacji
        act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
        act_constellation = setdiff(act_constellation, fail_sats); % usuniêcie 220
        nsat = length(act_constellation);
        
        for k=1:nsat
            
           prn_num = act_constellation(k);
           prn_idx = prn_num - sysPrefix;
           % interpolacja wsp satelity na epoke
           i_sp3 = find((eph(:,2)==prn_num));
           X_int = eph(i_sp3,1);
           Y_int = [eph(i_sp3,3) eph(i_sp3,4) eph(i_sp3,5)];          
           Xs = lagrange(X_int, Y_int, epochs(i), 10);
           % pêtla do obliczenia tau
           for s=1:2
                if (s > 1)
                   tau = (geom(k))/c;
                end               
                Xs  = e_r_corr(tau, Xs');            
                [azymut(k), wys(k), geom(k)] = topocent(X2, Xs-X2);              
                Xs = Xs';
           end
           % wsp sat po poprawce zegara odbiornika i tau
           Xs = lagrange(X_int, Y_int, epochs(i)+dtrec2(i,2)-tau, 10);
           % znowu obrót
           Xs  = e_r_corr(tau, Xs');
           %ostateczne wartoœci azymutu, wysokoœci, odleglosci i tropo
           [azymut1(k), wys1(k), geom1(k)] = topocent(X1, Xs-X1);
           [azymut2(k), wys2(k), geom2(k)] = topocent(X2, Xs-X2);
           dtropo1(k)=tropo(wys1(k), H1);
           dtropo2(k)=tropo(wys2(k), H2);
           
           obs1(k) = obsMatrix1(i_epoch, prn_idx, i_obs);
           obs2(k) = obsMatrix2(i_epoch, prn_idx, i_obs);
           % wyci¹gniêcie tylko czêœci wspólnej dla obu odbiorników
           if isnan(obs1(k))
               obs1(k)=[];
               obs2(k) = [];
               del_1(k) = k;
           elseif isnan(obs2(k))
               obs2(k) = [];
               obs1(k) = [];
               del_2(k) = k;
           end
           
           % satelita referencyjny 
           if k==nsat
               sat_ref1 = find_ref_sat(wys1);
               sat_ref2 = find_ref_sat(wys2);
               %sprawdzenie czy obserwuje wszystkie sygna³y
               if ismember(sat_ref1,del_1)
                   wys1_1(del_1) = [];
                   max_1 = max(wys1_1);
                   sat_ref1 = find(wys1==max_1);
               end
               if ismember(sat_ref2,del_2)
                   wys2_1(del_2) = [];
                   max_2 = max(wys2_1);
                   sat_ref2 = find(wys2==max_2);
               end
               % ustawienie na pierwszym miejscu w obs
               obs1([1 sat_ref1]) = obs1([sat_ref1 1]);
               obs2([1 sat_ref2]) = obs2([sat_ref2 1]);           
           end        
                      
        end
        
        % zapis obserwacji do odpowiedniej tablicy wed³ug typu
        if string(obs_types(j)) == "C1C"
            C1a = obs1; 
            C1b = obs2;            
        elseif string(obs_types(j)) == "L1C"
            L1a = obs1;
            L1b = obs2;     
        elseif string(obs_types(j)) == "C6C"
            C6a = obs1;
            C6b = obs2;       
        elseif string(obs_types(j)) == "L6C"
            L6a = obs1;
            L6b = obs2;      
        elseif string(obs_types(j)) == "C5Q"
            C5a = obs1;
            C5b = obs2;        
        elseif string(obs_types(j)) == "L5Q"
            L5a = obs1;
            L5b = obs2;        
        elseif string(obs_types(j)) == "C7Q"
           C7a = obs1;
           C7b = obs2;       
        elseif string(obs_types(j)) == "L7Q"
            L7a = obs1;
            L7b = obs2;        
        elseif string(obs_types(j)) == "C8Q"
            C8a = obs1;
            C8b = obs2;        
        elseif string(obs_types(j)) == "L8Q"
            L8a = obs1;
            L8b = obs2;
        end
        
        % czyszczenie przed nast iteracj¹
        obs1 = [];
        obs2 = [];
        
        
    end
    
end
