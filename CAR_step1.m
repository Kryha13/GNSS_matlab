function [Ncom1] = CAR_step1(obs_types, time_interval)
const;
load('dane_do_test.mat')
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
dtrec_idx = 1;
com1 = [0 0 1 -1];
com2 = [0 1 -1 0];
com3 = [1 -1 0 0];
com4 = [0 0 0 1];
D = [com1; com2; com3; com4];
phase_freq = [fE1 fE6 fE5 fE7];
% C = D^-1;

    for i=1:length(epochs)
        
        %% STEP 1
        
        [L, ~, DU, CL, obs_qty] = create_obs_matrices(epochs(i), obs_types, dtrec_idx);

        B = B_matrix_CAR(1, obs_qty);

        B = [DU B];
        xN = (B'*CL^-1*B)^-1 * B'*CL^-1*L;

        x_float(:,i) = xN(1:3);
        N_float = xN(4:3+str2double(obs_qty(1,2)));

        Cv = (B'*CL*B)^-1; 
        Cx = Cv(1:3,1:3);
        Cn = Cv(4:3+str2double(obs_qty(1,2)),4:3+str2double(obs_qty(1,2)));
    %     CxN = Cv(1:3, 4:length(xN));
    %     CNx = Cv(4:length(xN),1:3);

        [N_fixed,sqnorm,~,~,~,~,~]=LAMBDA(N_float, Cn);  %
        ratio(i) = sqnorm(:,2)./sqnorm(:,1);
        N_fixed1 = N_fixed(:,1); % pierwsza kolumna z Nfixed

        Ncom1 = N_fixed1;

        %% STEP 2
%            ["C1C" "L1C" "C6C" "L6C" "C5Q" "L5Q" "C7Q" "L7Q"];
        
        [L4, ~, DU34, CL34, obs_qty34] = create_obs_matrices(epochs(i), ["L7Q"], dtrec_idx);
        [L3, ~, ~, ~, ~] = create_obs_matrices(epochs(i), ["L5Q"], dtrec_idx);
        
        L4 = L4 * phase_freq(4)/c;
        L3 = L3 * phase_freq(3)/c;        
        [lam_com1] = get_com_lambda(1);
        Lcom1 = (L3 - L4)*lam_com1;
        
        [L, ~, DU, CL, obs_qty] = create_obs_matrices(epochs(i), ["L1C" "L6C" "L5Q"], dtrec_idx);
        
        obs_qty = [obs_qty34; obs_qty];
        
        dtrec_idx = dtrec_idx + 1;
    end
    
end
