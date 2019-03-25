const;
% time_interval = 100;
X1 = [3655333.1790; 1403901.2258; 5018038.3565]; %SEPT
X2 = [3655336.8160; 1403899.0382; 5018036.4527]; % SEP2
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
act_time1 = epochs + dtrec1(:,2);
act_time2 = epochs + dtrec2(:,2);
fodw = 298.2572221;
tau = 0.07; % przybli�ony 
% do stochastycznego
ac = 0.75; 
sigma = 0.002;
covL = 0.003;
covP = 0.35;

obs_types = ["C1C" "L1C" "C6C" "L6C" "C5Q" "L5Q" "C7Q" "L7Q" "C8Q" "L8Q"]; % wyb�r obserwacji 
phase_freq = [fE1 fE6 fE5 fE7 fE8];
system = 'E';
[sysPrefix, sysName] = sysGNSS(system);
fail_sats = [220 222]; % odbierany by� sygna� z 220 w obs a w efemerydach go nie ma nawet

[B1, L1, H1] = togeod(a, fodw, X1(1), X1(2), X1(3));
[B2, L2, H2] = togeod(a, fodw, X2(1), X2(2), X2(3));

% wybra� tylko te kt�re s� w obs_types
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
info = {};

for i=1:length(epochs)
    
    i_epoch = find(cell2mat(obsTable1(:,1))==epochs(i)); % indeks wiersza w obsMatrix
    
    for j=1:length(obs_types)
        
        i_obs = find(string(obsType1(:))==obs_types(j)); % indeks warstwy w obsMatrix
        % konstelacje dla obu odbiornik�w w danej epoce
        act_constellation_1 = cell2mat(obsTable1(i_epoch, i_obs+2)); 
        act_constellation_2 = cell2mat(obsTable2(i_epoch, i_obs+2));
        % cz�� wsp�lna konstelacji
        act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
        act_constellation = setdiff(act_constellation, fail_sats); % usuni�cie 220
        nsat = length(act_constellation);
        
        for k=1:nsat
            
           prn_num = act_constellation(k);
           prn_idx = prn_num - sysPrefix;
           % interpolacja wsp satelity na epoke
           i_sp3 = find((eph(:,2)==prn_num));
           X_int = eph(i_sp3,1);
           Y_int = [eph(i_sp3,3) eph(i_sp3,4) eph(i_sp3,5)];          
           Xs = lagrange(X_int, Y_int, epochs(i), 10);
           % p�tla do obliczenia tau
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
           % znowu obr�t
           Xs  = e_r_corr(tau, Xs');
           %ostateczne warto�ci azymutu, wysoko�ci, odleglosci i tropo
           [azymut1(k), wys1(k), geom1(k)] = topocent(X1, Xs-X1);
           [azymut2(k), wys2(k), geom2(k)] = topocent(X2, Xs-X2);
           dtropo1(k)=tropo(wys1(k), H1);
           dtropo2(k)=tropo(wys2(k), H2);
           u(k,:) = (Xs - X2)/geom2(k); % wersory u
           C(k) = (1+ac/sind(wys2(k)))^2 * sigma^2; % do stochastycznego
           
           obs1(k) = obsMatrix1(i_epoch, prn_idx, i_obs);
           obs2(k) = obsMatrix2(i_epoch, prn_idx, i_obs);
           % wyci�gni�cie tylko cz�ci wsp�lnej dla obu odbiornik�w
           if isnan(obs1(k))
               obs1(k)=[];
               obs2(k) = [];
               geom1(k) = [];
               geom2(k) = [];
               dtropo1(k) = [];
               dtropo2(k) = [];
               wys1(k) = [];
               wys2(k) = [];
               u(k) = [];
               C(k) = [];
               del_1(k) = k; %pomocnicze do sat ref
           elseif isnan(obs2(k))
               obs2(k) = [];
               obs1(k) = [];
               geom1(k) = [];
               geom2(k) = [];
               dtropo1(k) = [];
               dtropo2(k) = [];
               wys1(k) = [];
               wys2(k) = [];
               u(k) = [];
               C(k) = [];
               del_2(k) = k;
           end
           
           % satelita referencyjny 
           if k==nsat
               sat_ref1 = find_ref_sat(wys1);
               sat_ref2 = find_ref_sat(wys2);
               %sprawdzenie czy obserwuje wszystkie sygna�y
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
               % ustawienie na pierwszym miejscu w obserwacjach
               obs1([1 sat_ref1]) = obs1([sat_ref1 1]);
               obs2([1 sat_ref2]) = obs2([sat_ref2 1]);
               
               dtropo1([1 sat_ref1]) = dtropo1([sat_ref1 1]);
               dtropo2([1 sat_ref2]) = dtropo2([sat_ref2 1]);
               
               geom1([1 sat_ref1]) = geom1([sat_ref1 1]);
               geom2([1 sat_ref2]) = geom2([sat_ref2 1]);
               
               wys1([1 sat_ref1]) = wys1([sat_ref1 1]);
               wys2([1 sat_ref2]) = wys2([sat_ref2 1]);
               
               u([1 sat_ref2],:) = u([sat_ref2 1],:);
               C([1 sat_ref2]) = C([sat_ref2 1]);
               
           end        
                      
        end
        
        % zapis obserwacji do odpowiedniej tablicy wed�ug typu
        if string(obs_types(j)) == "C1C"
            C1a = obs1; 
            C1b = obs2;
            dtropo1C1 = dtropo1;
            dtropo2C1 = dtropo2;
            geomC1 = (geom2)';
            uC1 = u;
            CC1 = C*100;
        elseif string(obs_types(j)) == "L1C"
            L1a = obs1;
            L1b = obs2;
            dtropo1L1 = dtropo1;
            dtropo2L1 = dtropo2;
            geomL1 = (geom2)';
            uL1 = u;
            CL1 = C;
        elseif string(obs_types(j)) == "C6C"
            C6a = obs1;
            C6b = obs2;
            dtropo1C6 = dtropo1;
            dtropo2C6 = dtropo2;
            geomC6 = (geom2)';
            uC6 = u;
            CC6 = C*100;
        elseif string(obs_types(j)) == "L6C"
            L6a = obs1;
            L6b = obs2;
            dtropo1L6 = dtropo1;
            dtropo2L6 = dtropo2;
            geomL6 = (geom2)';
            uL6 = u;
            CL6 = phase_freq(1)/phase_freq(2) * C;
        elseif string(obs_types(j)) == "C5Q"
            C5a = obs1;
            C5b = obs2;
            dtropo1C5 = dtropo1;
            dtropo2C5 = dtropo2;
            geomC5 = (geom2)';
            uC5 = u;
            CC5 = C*100;
        elseif string(obs_types(j)) == "L5Q"
            L5a = obs1;
            L5b = obs2;
            dtropo1L5 = dtropo1;
            dtropo2L5 = dtropo2;
            geomL5 = (geom2)';
            uL5 = u;
            CL5 = phase_freq(1)/phase_freq(3) * C;
        elseif string(obs_types(j)) == "C7Q"
            C7a = obs1;
            C7b = obs2;
            dtropo1C7 = dtropo1;
            dtropo2C7 = dtropo2;
            geomC7 = (geom2)';
            uC7 = u;
            CC7 = C*100;
        elseif string(obs_types(j)) == "L7Q"
            L7a = obs1;
            L7b = obs2;
            dtropo1L7 = dtropo1;
            dtropo2L7 = dtropo2;
            geomL7 = (geom2)';
            uL7 = u;
            CL7 = phase_freq(1)/phase_freq(4) * C;
        elseif string(obs_types(j)) == "C8Q"
            C8a = obs1;
            C8b = obs2;
            dtropo1C8 = dtropo1;
            dtropo2C8 = dtropo2;
            geomC8 = (geom2)';
            uC8 = u;
            CC8 = C*100;
        elseif string(obs_types(j)) == "L8Q"
            L8a = obs1;
            L8b = obs2;
            dtropo1L8 = dtropo1;
            dtropo2L8 = dtropo2;
            geomL8 = (geom2)';
            uL8 = u;
            CL8 = phase_freq(1)/phase_freq(5) * C;
        end
        
        % czyszczenie przed nast iteracj�
        obs1 = [];
        obs2 = [];
        wys1 = [];
        wys2 = [];
        geom1 = [];
        geom2 = [];
        dtropo1 = [];
        dtropo2 = [];
        u = [];
        C = [];
        
    end
    % zapis podstawowych informacji dla epok
    info(1, :, i) = {epochs(i) sat_ref1 sat_ref2 act_constellation};
    
    % pojedyncze r�nice i macierze tworz�ce podw�jnych
    single_diff_C1 = (C1b - C1a)';
    single_diff_tropoC1 = (dtropo2C1 - dtropo1C1)';
    d_C1 = design_matrix(length(single_diff_C1));
    
    single_diff_L1 = (L1b - L1a)';
    single_diff_tropoL1 = (dtropo2L1 - dtropo1L1)';
    d_L1 = design_matrix(length(single_diff_L1));
    
    single_diff_C6 = (C6b - C6a)';
    single_diff_tropoC6 = (dtropo2C6 - dtropo1C6)';
    d_C6 = design_matrix(length(single_diff_C6));
    
    single_diff_L6 = (L6b - L6a)';
    single_diff_tropoL6 = (dtropo2L6 - dtropo1L6)';
    d_L6 = design_matrix(length(single_diff_L6));
    
    single_diff_C5 = (C5b - C5a)';
    single_diff_tropoC5 = (dtropo2C5 - dtropo1C5)';
    d_C5 = design_matrix(length(single_diff_C5));
    
    single_diff_L5 = (L5b - L5a)';
    single_diff_tropoL5 = (dtropo2L5 - dtropo1L5)';
    d_L5 = design_matrix(length(single_diff_L5));
    
    single_diff_C7 = (C7b - C7a)';
    single_diff_tropoC7 = (dtropo2C7 - dtropo1C7)';
    d_C7 = design_matrix(length(single_diff_C7));
    
    single_diff_L7 = (L7b - L7a)';
    single_diff_tropoL7 = (dtropo2L7 - dtropo1L7)';
    d_L7 = design_matrix(length(single_diff_L7));
    
    single_diff_C8 = (C8b - C8a)';
    single_diff_tropoC8 = (dtropo2C8 - dtropo1C8)';
    d_C8 = design_matrix(length(single_diff_C8));
    
    single_diff_L8 = (L8b - L8a)';
    single_diff_tropoL8 = (dtropo2L8 - dtropo1L8)';
    d_L8 = design_matrix(length(single_diff_L8));
    
    %podw�jne r�nice, wersory u, macierze B
    double_diff_C1 = d_C1*single_diff_C1;
    double_diff_tropoC1 = d_C1*single_diff_tropoC1;
    du_C1 = d_C1 * uC1;
    geom_C1 = d_C1 * geomC1; 
    
    double_diff_L1 = d_L1*single_diff_L1;
    double_diff_tropoL1 = d_L1*single_diff_tropoL1;
    du_L1 = d_L1 * uL1;
    geom_L1 = d_L1 * geomL1;
    B_L1 = eye(length(double_diff_L1))*phase_freq(1);
    
    double_diff_C6 = d_C6*single_diff_C6;
    double_diff_tropoC6 = d_C6*single_diff_tropoC6;
    du_C6 = d_C6 * uC6;
    geom_C6 = d_C6 * geomC6;
    
    double_diff_L6 = d_L6*single_diff_L6;
    double_diff_tropoL6 = d_L6*single_diff_tropoL6;
    du_L6 = d_L6 * uL6;
    geom_L6 = d_L6 * geomL6;
    B_L6 = eye(length(double_diff_L6))*phase_freq(2);
    
    double_diff_C5 = d_C5*single_diff_C5;
    double_diff_tropoC5 = d_C5*single_diff_tropoC5;
    du_C5 = d_C5 * uC5;
    geom_C5 = d_C5 * geomC5;
    
    double_diff_L5 = d_L5*single_diff_L5;
    double_diff_tropoL5 = d_L5*single_diff_tropoL5;
    du_L5 = d_L5 * uL5;
    geom_L5 = d_L5 * geomL5;
    B_L5 = eye(length(double_diff_L5))*phase_freq(3);
    
    double_diff_C7 = d_C7*single_diff_C7;
    double_diff_tropoC7 = d_C7*single_diff_tropoC7;
    du_C7 = d_C7 * uC7;
    geom_C7 = d_C7 * geomC7;
    
    double_diff_L7 = d_L7*single_diff_L7;
    double_diff_tropoL7 = d_L7*single_diff_tropoL7;
    du_L7 = d_L7 * uL7;
    geom_L7 = d_L7 * geomL7;
    B_L7 = eye(length(double_diff_L7))*phase_freq(4);
    
    double_diff_C8 = d_C8*single_diff_C8;
    double_diff_tropoC8 = d_C8*single_diff_tropoC8;
    du_C8 = d_C8 * uC8;
    geom_C8 = d_C8 * geomC8;
    
    double_diff_L8 = d_L8*single_diff_L8;
    double_diff_tropoL8 = d_L8*single_diff_tropoL8;
    du_L8 = d_L8 * uL8;
    geom_L8 = d_L8 * geomL8;
    B_L8 = eye(length(double_diff_L8))*phase_freq(5);
    
    %% zestawienie macierzy w modelu deterministycznym
    
    %wektor podw�jnych r�nic obs
    L = [double_diff_L1;double_diff_L6;double_diff_L5; double_diff_L7; double_diff_L8; ...
        double_diff_C1;double_diff_C6;double_diff_C5;double_diff_C7;double_diff_C8];
    T = [double_diff_tropoL1;double_diff_tropoL6;double_diff_tropoL5;double_diff_tropoL7;double_diff_tropoL8;...
        double_diff_tropoC1;double_diff_tropoC6;double_diff_tropoC5;double_diff_tropoC7;double_diff_tropoC8];
    % ?????
    R = [geom_L1;geom_L6;geom_L5;geom_L7;geom_L8;...
        geom_C1;geom_C6;geom_C5;geom_C7;geom_C8];
    
    L = L - R - T;
    
    % macierz wersor�w 
    DU = [du_L1;du_L6;du_L5;du_L7;du_L8;...
        du_C1;du_C6;du_C5;du_C7;du_C8];
    
    %macierz B
    B1 = blkdiag(B_L1, B_L6, B_L5, B_L7, B_L8);
    B0 = zeros(size(B1));
    B = [B1; B0];
    
    %% model stochastyczny
    
    % C0
    
    C0 = blkdiag(diag(CL1),diag(CL6),diag(CL5),diag(CL7),diag(CL8), ...
       diag(CC1), diag(CC6), diag(CC5), diag(CC7), diag(CC8));
    
    D = blkdiag(d_L1,d_L6,d_L5,d_L7,d_L8,d_C1,d_C6,d_C5,d_C7,d_C8);
    
    CL = 2 * D * C0 * D';
    
  
end
   






