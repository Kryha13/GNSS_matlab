function [com_lam] = get_com_lambda(combination)
    
const;
phase_freq = [fE1 fE6 fE5 fE7];
lam_up = c/fE1 * c/fE6 * c/fE5 * c/fE7;
lam_down = [c/fE6*c/fE5*c/fE7 c/fE1*c/fE5*c/fE7 c/fE1*c/fE6*c/fE7 c/fE1*c/fE6*c/fE5];
com1 = [0; 0; 1; -1];
com2 = [0; 1; -1; 0];
com3 = [1; -1; 0; 0];
com4 = [0; 0; 0; 1];

    switch combination       
        case 1
            com_lam = abs((lam_up/sum(lam_down*com1)));
        case 2
            com_lam = abs((lam_up/sum(lam_down*com2)));
        case 3 
            com_lam = abs((lam_up/sum(lam_down*com3)));
    end
end

