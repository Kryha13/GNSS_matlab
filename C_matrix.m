function [C] = C_matrix(Az_EL_R1, Az_EL_R2, time_interval, phase_frequencies, data)
load(data)

a = 0.75; 
sigma = 0.002;
covL = 0.003;
covP = 0.35;
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
n_obs = length(phase_frequencies)*2; 
n_sup = n_obs-1;
C = {};

for i=1:length(epochs)

    act_constellation_1 = squeeze(Az_EL_R1(:,2,i)); 
    act_constellation_2 = squeeze(Az_EL_R2(:,2,i));
    act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable');
    low_elev = find(Az_EL_R1(:, 4, i)<10);
    low_elev_satts = Az_EL_R1(low_elev, 2, i);
    act_constellation = setdiff(act_constellation, low_elev_satts);
    
    epoch_elev = Az_EL_R1(:,4,i);
    epoch_elev1 = Az_EL_R1(low_elev,4,i);
    epoch_elev = setdiff(epoch_elev, epoch_elev1);
    
    nsat=length(act_constellation);
    [sat_ref, sat_idx] = find_ref_sat(Az_EL_R1,i);
    sat_ref = find(act_constellation==sat_ref); %indeks w akt konstelacji
    
    for j=1:nsat
        
       CL1(j) = (1+a/sind(epoch_elev(j)))^2 * sigma^2;
       CL2(j) = phase_frequencies(1)/phase_frequencies(2) * CL1(j);
       CL3(j) = phase_frequencies(1)/phase_frequencies(3) * CL1(j);
       CL4(j) = phase_frequencies(1)/phase_frequencies(4) * CL1(j);
       CL5(j) = phase_frequencies(1)/phase_frequencies(5) * CL1(j);
       
       % satelita referencyjny jako pierwszy
        if j==nsat
            CL1([1 sat_ref]) = CL1([sat_ref 1]);
            CL2([1 sat_ref]) = CL2([sat_ref 1]);
            CL3([1 sat_ref]) = CL3([sat_ref 1]);
            CL4([1 sat_ref]) = CL4([sat_ref 1]);
            CL5([1 sat_ref]) = CL5([sat_ref 1]);
        end

       covL1L2(j) = phase_frequencies(1)/phase_frequencies(2) * covL^2;
       covL2L3(j) = phase_frequencies(2)/phase_frequencies(3) *covL^2;
       covL3L4(j) = phase_frequencies(3)/phase_frequencies(4) *covL^2;
       covL4L5(j) = phase_frequencies(4)/phase_frequencies(5) *covL^2; 
       
       covPk(j) = covP^2;
          
    end
   
   n = length(act_constellation);
   zero = zeros(n);
   Asup = [diag(covL1L2) diag(covL2L3) diag(covL3L4) diag(covL4L5) zero ...
           diag(covPk) diag(covPk) diag(covPk) diag(covPk)];
   Asup = reshape(Asup, [n n n_sup]); 
   Asub = Asup;
   Amd = [diag(CL1) diag(CL2) diag(CL3) diag(CL4) diag(CL5) ...
       diag(CL1)*100 diag(CL2)*100 diag(CL3)*100 diag(CL4)*100 diag(CL5)*100];
   Amd = reshape(Amd, [n n n_obs]);
   
   C{i} = full(blktridiag(Amd, Asub, Asup));
    
end

end

