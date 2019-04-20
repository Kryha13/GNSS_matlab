function [L, B, DU, CL] = create_obs_matrices(epoch, obs_types, dtrec_idx)
const;
load('dane_do_test.mat')
X1 = X(1,:)'; %SEPT
X2 = X(2,:)'; % SEP2
epochs = epoch;
fodw = 298.2572221;
tau1 = 0.07;% przybli¿ony
tau2 = 0.07;
% do stochastycznego
ac = 0.75; 
sigma = 0.002;

phase_freq = [fE1 fE6 fE5 fE7 fE8];
system = 'E';
[sysPrefix, ~] = sysGNSS(system);
fail_sats = [20 22]; 
[~, ~, H1] = togeod(a, fodw, X1(1), X1(2), X1(3));
[~, ~, H2] = togeod(a, fodw, X2(1), X2(2), X2(3));



for i=1:length(epochs)
    L_L_all = [];
    L_c_all = [];
    B = [];
    DU_L = [];
    DU_c = [];
    CL_L = [];
    CL_c = [];
    
    i_epoch = find(cell2mat(obsTable1(:,1))==epochs(i)); % indeks wiersza w obsMatrix
    
    for j=1:length(obs_types)
        
        i_obs = find(string(obsType1(:))==obs_types(j)); % indeks warstwy w obsMatrix
        % konstelacje dla obu odbiorników w danej epoce
        act_constellation_1 = find(~isnan(obsMatrix1(i_epoch,:,i_obs))); 
        act_constellation_2 = find(~isnan(obsMatrix2(i_epoch,:,i_obs)));
        % czêœæ wspólna konstelacji
        act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
        act_constellation = setdiff(act_constellation, fail_sats); % usuniêcie fail
        nsat = length(act_constellation);
        
        for k=1:nsat
            
           prn_idx = act_constellation(k);
           prn_num = act_constellation(k) + sysPrefix;
           % interpolacja wsp satelity na epoke
           i_sp3 = find((eph(:,2)==prn_num));
           X_int = eph(i_sp3,1);
           Y_int = [eph(i_sp3,3) eph(i_sp3,4) eph(i_sp3,5)];          
           Xs = lagrange(X_int, Y_int, epochs(i), 10);
           % pêtla do obliczenia tau
           for s=1:2
                if (s > 1)
                   tau1 = (geo1(k))/c;
                   tau2 = (geo2(k))/c;
                end               
                Xs1  = e_r_corr(tau1, Xs');
                Xs2  = e_r_corr(tau2, Xs');
                [az1(k), w1(k), geo1(k)] = topocent(X1, Xs1-X1);
                [az2(k), w2(k), geo2(k)] = topocent(X2, Xs2-X2);              
                Xs1 = Xs1';
                Xs2 = Xs2';
           end
           Xs1 = [];
           Xs2 = [];
           % wsp sat po poprawce zegara odbiornika i tau
           Xs1 = lagrange(X_int, Y_int, epochs(i)-dtrec1(dtrec_idx,2)-tau1, 10);
           Xs2 = lagrange(X_int, Y_int, epochs(i)-dtrec2(dtrec_idx,2)-tau2, 10);
           % znowu obrót
           Xs1  = e_r_corr(tau1, Xs1');
           Xs2 = e_r_corr(tau2, Xs2');
           %ostateczne wartoœci azymutu, wysokoœci, odleglosci i tropo
           [azymut1(k), wys1(k), geom1(k)] = topocent(X1, Xs1-X1);
           [azymut2(k), wys2(k), geom2(k)] = topocent(X2, Xs2-X2);
           dtropo1(k)=tropo(wys1(k), H1);
           dtropo2(k)=tropo(wys2(k), H2);
           u(k,:) = (Xs2 - X2)/geom2(k); % wersory u
           C(k) = (1+ac/sind(wys2(k)))^2 * sigma^2; % do stochastycznego
           
           obs1(k) = obsMatrix1(i_epoch, prn_idx, i_obs);
           obs2(k) = obsMatrix2(i_epoch, prn_idx, i_obs);
           
           % satelita referencyjny 
           if k==nsat
               sat_ref1 = find_ref_sat(wys1);
               sat_ref2 = find_ref_sat(wys2);
               %sprawdzenie czy obserwuje wszystkie sygna³y
               % ustawienie na pierwszym miejscu w obserwacjach
               while 1                 
                   if any(isnan(obsMatrix1(i_epoch, sat_ref1,:)))
                       wys11 = wys1;
                       wys11(sat_ref1) = [];
                       sat_ref1 = find_ref_sat(wys11);
                       sat_ref1 = find(wys1==wys11(sat_ref1));                       
                   elseif any(isnan(obsMatrix2(i_epoch, sat_ref2,:)))
                       wys22 = wys2;
                       wys22(sat_ref2) = [];
                       sat_ref2 = find_ref_sat(wys22);
                       sat_ref2 = find(wys2==wys22(sat_ref2));
                   end                  
                   if all(~isnan(obsMatrix1(i_epoch, sat_ref1,:))) && all(~isnan(obsMatrix2(i_epoch, sat_ref2,:)))
                       break
                   end
                end
               
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
                      
        % zapis obserwacji do odpowiedniej tablicy wed³ug typu
        if string(obs_types(j)) == "C1C" || string(obs_types(j)) == "L1C"
            % CODE
            if string(obs_types(j)) == "C1C"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_c1 = C*2500; %%%%%
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_c = design_matrix(length(single_diff)); %%%%%%
                double_diff = d_c*single_diff;
                double_diff_tropo = d_c*single_diff_tropo;
                double_diff_geom = d_c*single_diff_geom;
                du_c1 = d_c * u; 

                L_c = [double_diff];
                T_c = [double_diff_tropo];
                R_c = [double_diff_geom];
                L_c1 = L_c - R_c - T_c;
                
                C0 = diag(C_c1);

                Cc_1 = 2 * d_c * C0 * d_c';
                
                L_c_all = [L_c_all; L_c1];
                DU_c = [DU_c; du_c1];
                CL_c = blkdiag(CL_c, Cc_1);

            elseif string(obs_types(j)) == "L1C"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_L1 = C;
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_L = design_matrix(length(single_diff));
                double_diff = d_L*single_diff*(c/phase_freq(1));
                double_diff_tropo = d_L*single_diff_tropo;
                double_diff_geom = d_L*single_diff_geom;
                du_L1 = d_L * u; 

                L_L = [double_diff];
                T_L = [double_diff_tropo];
                % ?????
                R_L = [double_diff_geom];
                %macierz B
                B_L = eye(length(double_diff))*(c/phase_freq(1));
                const_L = act_constellation;
                
                L_1 = L_L - R_L - T_L;
                B1 = blkdiag(B_L);
                B_1 = [B1];

                C0 = diag(C_L1);

                CL_1 = 2 * d_L * C0 * d_L';
                
                L_L_all = [L_L_all; L_1];
                B = blkdiag(B, B_1);
                DU_L = [DU_L; du_L1];
                CL_L = blkdiag(CL_L, CL_1);          
            end


        elseif string(obs_types(j)) == "C6C" || string(obs_types(j)) == "L6C"
            % CODE
            if string(obs_types(j)) == "C6C"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_c2 = C*2500; %%%%%
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_c = design_matrix(length(single_diff)); %%%%%%
                double_diff = d_c*single_diff;
                double_diff_tropo = d_c*single_diff_tropo;
                double_diff_geom = d_c*single_diff_geom;
                du_c2 = d_c * u; 

                L_c = [double_diff];
                T_c = [double_diff_tropo];
                R_c = [double_diff_geom];
                L_c2 = L_c - R_c - T_c;
                
                C0 = diag(C_c2);

                Cc_2 = 2 * d_c * C0 * d_c';
                
                L_c_all = [L_c_all; L_c2];
                DU_c = [DU_c; du_c2];
                CL_c = blkdiag(CL_c, Cc_2);
                               
            elseif string(obs_types(j)) == "L6C"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_L2 = C;
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_L = design_matrix(length(single_diff));
                double_diff = d_L*single_diff*(c/phase_freq(2));
                double_diff_tropo = d_L*single_diff_tropo;
                double_diff_geom = d_L*single_diff_geom;
                du_L2 = d_L * u; 

                L_L = [double_diff];
                T_L = [double_diff_tropo];
                % ?????
                R_L = [double_diff_geom];
                %macierz B
                B_L = eye(length(double_diff))*(c/phase_freq(2));
                const_L = act_constellation;
                
                L_2 = L_L - R_L - T_L;
                B1 = blkdiag(B_L);
                B_2 = [B1];

                C0 = diag(C_L2);

                CL_2 = 2 * d_L * C0 * d_L';
                
                L_L_all = [L_L_all; L_2];
                B = blkdiag(B, B_2);
                DU_L = [DU_L; du_L2];
                CL_L = blkdiag(CL_L, CL_2); 
            end
            
        elseif string(obs_types(j)) == "C5Q" || string(obs_types(j)) == "L5Q"
            % CODE
            if string(obs_types(j)) == "C5Q"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_c3 = C*2500; %%%%%
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_c = design_matrix(length(single_diff)); %%%%%%
                double_diff = d_c*single_diff;
                double_diff_tropo = d_c*single_diff_tropo;
                double_diff_geom = d_c*single_diff_geom;
                du_c3 = d_c * u; 

                L_c = [double_diff];
                T_c = [double_diff_tropo];
                R_c = [double_diff_geom];
                L_c3 = L_c - R_c - T_c;
                
                C0 = diag(C_c3);

                Cc_3 = 2 * d_c * C0 * d_c';
                
                L_c_all = [L_c_all; L_c3];
                DU_c = [DU_c; du_c3];
                CL_c = blkdiag(CL_c, Cc_3);
                               
            elseif string(obs_types(j)) == "L5Q"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_L3 = C;
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_L = design_matrix(length(single_diff));
                double_diff = d_L*single_diff*(c/phase_freq(3));
                double_diff_tropo = d_L*single_diff_tropo;
                double_diff_geom = d_L*single_diff_geom;
                du_L3 = d_L * u; 

                L_L = [double_diff];
                T_L = [double_diff_tropo];
                % ?????
                R_L = [double_diff_geom];
                %macierz B
                B_L = eye(length(double_diff))*(c/phase_freq(3));
                const_L = act_constellation;
                
                L_3 = L_L - R_L - T_L;
                B1 = blkdiag(B_L);
                B_3 = [B1];

                C0 = diag(C_L3);

                CL_3 = 2 * d_L * C0 * d_L';
                
                L_L_all = [L_L_all; L_3];
                B = blkdiag(B, B_3);
                DU_L = [DU_L; du_L3];
                CL_L = blkdiag(CL_L, CL_3); 
                
            end
            
        elseif string(obs_types(j)) == "C7Q" || string(obs_types(j)) == "L7Q"
            % CODE
            if string(obs_types(j)) == "C7Q"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_c4 = C*2500; %%%%%
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_c = design_matrix(length(single_diff)); %%%%%%
                double_diff = d_c*single_diff;
                double_diff_tropo = d_c*single_diff_tropo;
                double_diff_geom = d_c*single_diff_geom;
                du_c4 = d_c * u; 

                L_c = [double_diff];
                T_c = [double_diff_tropo];
                R_c = [double_diff_geom];
                L_c4 = L_c - R_c - T_c;
                
                C0 = diag(C_c4);

                Cc_4 = 2 * d_c * C0 * d_c';
                
                L_c_all = [L_c_all; L_c4];
                DU_c = [DU_c; du_c4];
                CL_c = blkdiag(CL_c, Cc_4);
                               
            elseif string(obs_types(j)) == "L7Q"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_L4 = C;
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_L = design_matrix(length(single_diff));
                double_diff = d_L*single_diff*(c/phase_freq(4));
                double_diff_tropo = d_L*single_diff_tropo;
                double_diff_geom = d_L*single_diff_geom;
                du_L4 = d_L * u; 

                L_L = [double_diff];
                T_L = [double_diff_tropo];
                % ?????
                R_L = [double_diff_geom];
                %macierz B
                B_L = eye(length(double_diff))*(c/phase_freq(4));
                const_L = act_constellation;
                
                L_4 = L_L - R_L - T_L;
                B1 = blkdiag(B_L);
                B_4 = [B1];

                C0 = diag(C_L4);

                CL_4 = 2 * d_L * C0 * d_L';
                
                L_L_all = [L_L_all; L_4];
                B = blkdiag(B, B_4);
                DU_L = [DU_L; du_L4];
                CL_L = blkdiag(CL_L, CL_4); 
                
            end
            
        elseif string(obs_types(j)) == "C8Q" || string(obs_types(j)) == "L8Q"
            % CODE
            if string(obs_types(j)) == "C8Q"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_c5 = C*2500; %%%%%
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_c = design_matrix(length(single_diff)); %%%%%%
                double_diff = d_c*single_diff;
                double_diff_tropo = d_c*single_diff_tropo;
                double_diff_geom = d_c*single_diff_geom;
                du_c5 = d_c * u; 

                L_c = [double_diff];
                T_c = [double_diff_tropo];
                R_c = [double_diff_geom];
                L_c5 = L_c - R_c - T_c;
                
                C0 = diag(C_c5);

                Cc_5 = 2 * d_c * C0 * d_c';
                
                L_c_all = [L_c_all; L_c5];
                DU_c = [DU_c; du_c5];
                CL_c = blkdiag(CL_c, Cc_5);
                               
            elseif string(obs_types(j)) == "L8Q"
                a = obs1; 
                b = obs2;
                dtropo1 = dtropo1;
                dtropo2 = dtropo2;
                geom1 = geom1;
                geom2 = geom2;
                u = u;
                C_L5 = C;
                single_diff = (b - a)';
                single_diff_tropo = (dtropo2 - dtropo1)';
                single_diff_geom = (geom2 - geom1)';
                d_L = design_matrix(length(single_diff));
                double_diff = d_L*single_diff*(c/phase_freq(5));
                double_diff_tropo = d_L*single_diff_tropo;
                double_diff_geom = d_L*single_diff_geom;
                du_L5 = d_L * u; 

                L_L = [double_diff];
                T_L = [double_diff_tropo];
                % ?????
                R_L = [double_diff_geom];
                %macierz B
                B_L = eye(length(double_diff))*(c/phase_freq(5));
                const_L = act_constellation;
                
                L_5 = L_L - R_L - T_L;
                B1 = blkdiag(B_L);
                B_5 = [B1];

                C0 = diag(C_L5);

                CL_5 = 2 * d_L * C0 * d_L';
                
                L_L_all = [L_L_all; L_5];
                B = blkdiag(B, B_5);
                DU_L = [DU_L; du_L5];
                CL_L = blkdiag(CL_L, CL_5); 
                
            end
        end
        
        % czyszczenie przed nast iteracj¹
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

end

L = [L_L_all; L_c_all];
B = blkdiag(B);
B0 = zeros(size(B));
B = [B; B0];
DU = [DU_L; DU_c];
CL = blkdiag(CL_L, CL_c);


end

