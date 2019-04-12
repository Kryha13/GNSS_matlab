function [design_mat] = design_matrix(nsat)

d =zeros(nsat-1,nsat);

idx = 2:1:nsat;
    for i=1:(nsat-1)
        d(i,1) = 1;
        d(i,idx(i)) = -1;      
    end  
    design_mat = d;
%     design_mat = kron(eye(obs_count),d); % dla wielu na raz
    
end

