function [design_mat] = design_matrix(obs_count, nsat)

d =zeros(nsat-1,nsat);
design_mat=[];
idx = 2:1:nsat;
    for i=1:(nsat-1)
        d(i,1) = 1;
        d(i,idx(i)) = -1;      
    end   
    design_mat = kron(eye(obs_count),d);
end

