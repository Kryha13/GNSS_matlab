function [X_code, DOP, dtrec_mat, Az_EL_R, dtrop, tau_ost] = codepos(Xapr, codes, time_interval, system, data)
%%
%%%%INPUT%%%%%
% Xapr - wspó³rzêdne przybli¿one z pliku rinex
% codes - kody wybrane do kombinacji liniowej [str1 str2]
% time_interval - interwa³ miêdzy pierwsz¹ a ostatni¹ epok¹ dla którego
% chcemy otrzymywaæ pozycjê
% system - identyfikator systemu nawigacyjnego
% data - opcjonalnie - wczytane dane z plików rinex, sp3, clk [str1 str2..]

%%%%OUTPUT%%%%
% X_code - wektory wspó³rzêdnych kodowych dla kolejnych epok
% DOP - wspó³czynniki DOP dla kolejnych epok
% dtrec_mat - poprawki zegara odbiornika dla kolejnych epok
% Az_EL_R - azymuty, elewacje i odlegloœci do satelitów; warstwy to kolejne
% epoki
%%
load(data); %opcjonalnie
const  % wczytanie sta³ych

time=gpsSecondsFirst:time_interval:gpsSecondsLast;  % interwa³ na jaki chcemy pozycje
[sysPrefix, sysName] = sysGNSS(system); % dla Galileo 
fodw = 298.2572221;
tau=0.07; %przybli¿ony czas propagacji sygna³u
dtrec=0;
X = [];
% constellation = {};
% Wspolczynniki do kombinacji liniowej
f1 = get_frequency(codes(1));
f2 = get_frequency(codes(2));
alfa13= f1^2/(f1^2-f2^2);
alfa23=-f2^2/(f1^2-f2^2);
%% indeksy wartw w obsMatrix dla u¿ywanych kodów
i_code1 = find(string(obsType(:))==codes(1));
i_code2 = find(string(obsType(:))==codes(2));


    for i=1:length(time) % pêtla po czasie na który wyznaczamy pozyccje
        X=Xapr';
        i_epoch = find(cell2mat(obsTable(:,1))==time(i)); % indeks wiersza epoki dla obsMatrix
        act_constellation_1 = cell2mat(obsTable(i_epoch, i_code1+2)); 
        act_constellation_2 = cell2mat(obsTable(i_epoch, i_code2+2));
        % lista satelitów obserwowanych w epoce dla obu kodów
        act_constellation = intersect(act_constellation_1, act_constellation_2, 'stable'); 
      
        for j=1:3

           [B, L, H] = togeod(a, fodw, X(1), X(2), X(3));
           nsat=length(act_constellation);
           isat=0;

            for k=1:nsat % pêtla po satelitach z danej epoki

                prn_num = act_constellation(k); % numer k-ty satelita w danej epoce
                prn_idx = prn_num - sysPrefix; % index dla kolumny w obsMatrix

                P12 =[obsMatrix(i_epoch, prn_idx, i_code1) obsMatrix(i_epoch, prn_idx, i_code2)];
                P3= P12(:,1)*alfa13 + P12(:,2)*alfa23;

                if (j > 1)
                   tau = (geom(k))/c;
                end
                
                if j==3
                    tau_ost(k,:,i) = [time(i) prn_num tau];
                end

                i_sp3 = find((eph(:,2)==prn_num));
                i_clk = find((sat_clk(:,2)==prn_num));

                X_int = eph(i_sp3,1);
                Y_int = [eph(i_sp3,3) eph(i_sp3,4) eph(i_sp3,5)];
                Xs = lagrange(X_int, Y_int, time(i), 10); % lagrange 10 stopnia - wsp sat

                X_int_clk = sat_clk(i_clk,1);
                Y_int_clk = sat_clk(i_clk,3);
                dtsat = lagrange(X_int_clk, Y_int_clk, time(i), 2); % poprawka zegara satelity

                Xs1 = lagrange(X_int, Y_int, time(i) + 1, 10);  %wsp w kolejnej epoce, do obliczenia predkosci satelity 

                Xsat  = e_r_corr(tau, Xs');
                Xsat1  = e_r_corr(tau, Xs1');
                V = (Xsat1 - Xsat) / 1;  % prêdkoœæ satelity

                [azymut(k), wys(k), geom(k)] = topocent(X, Xsat-X);

                Az_EL_R(k,:,i) = [time(i) prn_num azymut(k) wys(k) geom(k)]; % warstwy to kolejne epoki

                if wys(k)>=10  % maska 10
%                     constellation(k,i) = {act_constellation(k)};
                   isat=isat+1;
                   PACT(isat)=P3;

                    % poprawki do pseudoodleg³oœci
%                    drel(k)= -(2*Xsat'*V)/c;
                   dtropo(k)=tropo(wys(k), H);
                   dtrop(k,:,i) = [time(i) prn_num dtropo(k)];

                   PCOM(isat) = geom(k) - dtsat*c  + dtropo(k);
                   A(isat,1) = -(Xsat(1) - X(1))/geom(k);
                   A(isat,2) = -(Xsat(2) - X(2))/geom(k);
                   A(isat,3) = -(Xsat(3) - X(3))/geom(k);
                   A(isat,4) = 1;
                end
            end % Koniec petli po satelitach 
            
           b= PACT' - PCOM';
           x=inv(A'*A)*(A'*b);
           v=A*x-b;
           m0=sqrt(v'*v/(isat-4));
           dtrec = (dtrec+x(4))/c;
           dtrec_mat(i,:) = [time(i) dtrec]; % poprawka zegara odbiornika
           X_code(i,:) = [time(i) (X+x(1:3))'];  % wspolrzedne

        end % Koniec petli jj

    Q = inv(A'*A);
    diagonal=diag(Q);
    GDOP = sqrt(sum(diagonal));
    PDOP = sqrt(diagonal(1)+diagonal(2)+diagonal(3));
    TDOP = sqrt(diagonal(4));
    DOP(i,:) = [time(i) GDOP PDOP TDOP];  % DOPY

    end % koniec petli po czasie
end
			