function [B] = B_matrix_CAR(step_number, obs_qty)

const;
phase_freq = [fE1 fE6 fE5 fE7];
lam_up = c/fE1 * c/fE6 * c/fE5 * c/fE7;
lam_down = [c/fE6*c/fE5*c/fE7 c/fE1*c/fE5*c/fE7 c/fE1*c/fE6*c/fE7 c/fE1*c/fE6*c/fE5];
com1 = [0; 0; 1; -1];
com2 = [0; 1; -1; 0];
com3 = [1; -1; 0; 0];
com4 = [0; 0; 0; 1];
D = [com1; com2; com3; com4];
phase_obs_qty = obs_qty(1:length(obs_qty)/2, 2);

    switch step_number        
        case 1
            if phase_obs_qty(3) == phase_obs_qty(4) && phase_obs_qty(1) == phase_obs_qty(2)
               B1 = eye(str2double(phase_obs_qty(1)))*(c/phase_freq(1));
               B2 = eye(str2double(phase_obs_qty(2)))*(c/phase_freq(2));
               B3 = eye(str2double(phase_obs_qty(3)))*(c/phase_freq(3));
               B4 = eye(str2double(phase_obs_qty(4)))*(c/phase_freq(4));
               B_L = [B1 B1 B1 B1; B2 B2 zeros(size(B2)) B2; ...
                   B3 zeros(size(B3)) zeros(size(B3)) B3;...
                   zeros(size(B4)) zeros(size(B4)) zeros(size(B4)) B4];
               B0 = zeros(size(B_L));
               B = [B_L; B0];
            else
                B = [];
            end
            
        case 2            
           if phase_obs_qty(3) == phase_obs_qty(4) && phase_obs_qty(1) == phase_obs_qty(2)
               Bcom2 = eye(str2double(phase_obs_qty(3)))*((lam_up/sum(lam_down*com2)));
               Bcom3 = eye(str2double(phase_obs_qty(1)))*((lam_up/sum(lam_down*com3)));
               B4 = eye(str2double(phase_obs_qty(4)))*(c/phase_freq(4));
               B_L = [Bcom2 Bcom3 B4; Bcom2 zeros(size(Bcom2)) B4; ...
                   zeros(size(Bcom2)) zeros(size(Bcom2)) B4];
               B0 = zeros(str2double(phase_obs_qty(1)));
               B = [B_L; B0];
           else
               B = [];
           end
           
        case 3
            if phase_obs_qty(3) == phase_obs_qty(4)
               Bcom3 = eye(str2double(phase_obs_qty(1)))*((lam_up/sum(lam_down*com3)));
               B4 = eye(str2double(phase_obs_qty(4)))*(c/phase_freq(4));
               B_L = [Bcom3 B4; zeros(size(Bcom3)) B4; zeros(size(Bcom3)) B4];
               B0 = zeros(str2double(phase_obs_qty(1)));
               B = [B_L; B0];
            end
            
        case 4
            B4 = eye(str2double(phase_obs_qty(4)))*(c/phase_freq(4));
            B_L = [B4; B4; B4];
            B0 = zeros(size(B_L));
            B = [B_L; B0];
        
    end
end

