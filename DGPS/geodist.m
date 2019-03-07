% geometric distance -----------------------------------------------------------
function [r,e,el,dt]=geodist(t0,t,rr,maska,nav)

%t0: data
%t= (epoka obs-poprawka zegara (x(4))/C)
%rr: wektor pozycji x(1:3)
%nav: efemeryda dla danego satelity - mo¿e byæ dla kilku epok

C=299792458; OMGE=7.292115167E-5; r=0; rk=1;
while abs(r-rk)>1E-4
    [rs,dt]=satpos(t0,t-r/C,nav);
    % rs: pozycja sat na epokê t-r/C
    % dt: poprawka zegara satelity
    
    rrs=rr-Rz(OMGE*r/C)*rs; rk=r; r=norm(rrs);
end
e=rrs/r;
if norm(rr)>0, el=-asin(rr'*e/norm(rr)); else el=pi/2; end
if el*180/pi<maska, r=nan; end % elevation cutoff