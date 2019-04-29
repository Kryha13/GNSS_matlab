function [x_fixed, ratio] = calculate_position(obs_types, time_interval)

load('dane_do_test.mat')
epochs = gpsSecondsFirst:time_interval:gpsSecondsLast;
dtrec_idx = 1;

for i=1:length(epochs)
    
    [L, B, DU, CL, ~] = create_obs_matrices(epochs(i), obs_types, dtrec_idx);

    %% rozwi¹zanie modelu pozycjonowania - FLOAT
    ATCA = DU' * (CL)^-1 * DU;
    ATCB = DU' * (CL)^-1 * B;
    BTCA = B' * (CL)^-1 * DU;
    BTCB = B' * (CL)^-1 * B;

    ATCL = DU' * (CL)^-1 * L;
    BTCL = B' * (CL)^-1 * L;

    M1 = [ATCA ATCB; BTCA BTCB];
    M2 = [ATCL; BTCL];

    % x_float
    xN = M1^-1 * M2;
    x_float(:,i) = xN(1:3);
    N_float = xN(4:length(xN));

    % macierze wariancyjno-kowariancyjne
    Cv = M1^-1; 
    Cx = Cv(1:3,1:3);
    Cn = Cv(4:length(xN),4:length(xN));
    CxN = Cv(1:3, 4:length(xN));
    CNx = Cv(4:length(xN),1:3);

    Dx_fixed = Cx - CxN*Cn^-1 * CNx;
    % LAMBDA
    [N_fixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(N_float,Cn);  %

    N_fixed1 = N_fixed(:,1); % pierwsza kolumna z Nfixed

    % obliczenei x fixed
    x_fixed(:,i) = x_float(:,i) - CxN*Cn^-1 * (N_float - N_fixed1);
    % ratio 
    ratio(i) = sqnorm(:,2)./sqnorm(:,1);
    
    dtrec_idx = dtrec_idx + 1;
    
     %% zerowanie 
    L = [];
    B = [];
    DU =[];
    CL = [];
    N_float = [];
    xN = [];
    Cv = [];
%     x_float = [];
    CxN = [];
    Cn = [];
%     N_float = [];
    N_fixed1 = [];
    Dx_fixed = [];
    
    ATCA = [];
    ATCB = [];
    BTCA = [];
    BTCB = [];
    
    ATCL = [];
    BTCL = [];
    
    M1 = [];
    M2 = [];
end
end

