function [x_fixed, ratio] = calculate_position(L, DU, B, CL)

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
    x_float = xN(1:3);
    N_float = xN(4:length(xN));
    
    % macierze wariancyjno-kowariancyjne
    Cv = M1^-1; 
    Cx = Cv(1:3,1:3);
    Cn = Cv(4:length(xN),4:length(xN));
    CxN = Cv(1:3, 4:length(xN));
    CNx = Cv(4:length(xN),1:3);
    
    Dx_fixed = Cx - CxN*Cn^-1 * CNx;
    % LAMBDA
    [N_fixed,sqnorm(i,:),Ps(i),Qzhat,Z,nfixed(i),mu(i)]=LAMBDA(N_float,Cn);  %
    
    N_fixed1 = N_fixed(:,1); % pierwsza kolumna z Nfixed
    
    % obliczenei x fixed
    x_fixed(:,i) = x_float(:,i) - CxN*Cn^-1 * (N_float - N_fixed1);
    % ratio 
    ratio = sqnorm(:,2)./sqnorm(:,1);
    
end

