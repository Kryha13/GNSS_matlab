function [X, DOP, dtrec,   ] = codepos(f1, f2, obs1, obs2, act_constellation, Xapr, time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% obliczenie pozycji kodowej %%
% INPUT:
%%% 
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fodw = 298.2572221;
nObs = 1;

P12 =[obs1 obs2];
% Wspolczynniki do kombinacji liniowej
alfa13= f1^2/(f1^2-f2^2);
alfa23=-f2^2/(f1^2-f2^2);
P3= P12(:,1)*alfa13 + P12(:,2)*alfa23 ;  %2.546*Bp*10^-9*c %człon z Bp kodowe opoznienie sprzętowe

% czas obserwacji, w sekundach tygodnia GPS
tobs=time; % nalezy policzyc dla swojej epoki, tu dla pierwszej
epoch=tobs;
% przyblizony czas propagacji sygnalu
tau=0.07;
dtrec=0;

X=Xapr';

	for jj=1:3

	   [B, L, H] = togeod(a, fodw, X(1), X(2), X(3));

	   [nsat,ncol]=size(act_constellation);

	   isat=0;

		for ii=1:nsat % petla po wszystkich satelitach

		   if (jj > 1)
			   tau = (geom(ii))/c;
		   end

	% wczytanie do wektora 'Xs' wspolrzednych kolejnych satelitow na moment tobs-tau
	% w zmiennej dtsat znajduje sie poprawka do zeg. satelity
		   % [Xs, dtsat]   = wspsat(act_constellation(ii), eph, tobs+dtrec-tau);
		   % [Xs1, dtsat1]   = wspsat(act_constellation(ii), eph, tobs+dtrec-tau +0.1);
			i_sp3 = find((eph(:,2)==prn_num));
			i_clk = find((sat_clk(:,2)==prn_num)) 
			
			X_int = eph(i_sp3,1);
			Y_int = [eph(i_sp3,3) eph(i_sp3,4) eph(i_sp3,5)];
			Xs = lagrange(X_int, Y_int, time(i), 10); % lagrange 10 stopnia - wsp sat
			
			X_int_clk = sat_clk(i_clk,1);
			Y_int_clk = sat_clk(i_clk,3);
			dtsat = lagrange(X_int_clk, Y_int_clk, time{i), 2); % poprawka zegara satelity
			
			Xs1 = lagrange(X_int, Y_int, time(i) + 1, 10)  %wsp w kolejnej epoce, do obliczenia predkosci satelity 
			
			% moga byc potrzbne transpozy wsp ???
			
	% Korekcja wsp. satelity z powodu obrotu Ziemi, skorzystaÄ‡ z procedury e_r_corr(tau, Xs)
		   Xsat  = e_r_corr(tau, Xs);
		   Xsat1  = e_r_corr(tau, Xs1);
		   
	% !!!! PoliczyÄ‡ prÄ™dkoĹ›Ä‡ satelity
		   V = (Xsat1 - Xsat) / 1

	% azymut, wysokosc nad horyzontem oraz odleglosc do satelity
		   [azymut(ii), wys(ii), geom(ii)] = topocent(X, Xsat-X);
		   if wys(ii)>=10
			   isat=isat+1;
			   PACT(isat)=P3(ii);
	% !!!! Poprawki do pseudoodlegĹ‚oĹ›ci, ktĂłre naleĹĽy uwzglÄ™dniÄ‡:
			   drel(ii)= -(2*Xsat'*V)/c;
			   dtropo(ii)=tropo(wys(ii), H);
			   dtrec;
	%###################
			   PCOM(isat) = geom(ii) - dtsat*c - drel(ii) + dtropo(ii);
			   A(isat,1) = -(Xsat(1) - X(1))/geom(ii);
			   A(isat,2) = -(Xsat(2) - X(2))/geom(ii);
			   A(isat,3) = -(Xsat(3) - X(3))/geom(ii);
			   A(isat,4) = 1;
			end
	   end % Koniec petli ii=1:nsat
	   
	   b= PACT' - PCOM';
	   x=inv(A'*A)*(A'*b);
	   v=A*x-b;
	   m0=sqrt(v'*v/(isat-4));
	   dtrec = dtrec+x(4)/c;
	   X=X+x(1:3);
	   
	end % Koniec petli jj


%DOPy, 

Q = inv(A'*A);
diag=diag(Q);

GDOP = sqrt(sum(diag));
PDOP = sqrt(diag(1)+diag(2)+diag(3));
TDOP = sqrt(diag(4));
DOP = [GDOP PDOP TDOP];

end
