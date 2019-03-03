% navigation messeage to satellite position/clock bias -------------------------
function [rs,dt]=satpos(t0,t,nav)
if isempty(nav), rs=[nan;nan;nan]; dt=0; return, end
C=299792458; MU=3.986005E14; OMGE=7.292115167E-5;
tk=(t0-datenum(1980,1,6)-nav(:,28)*7)*86400+t-nav(:,18);
[tt,i]=min(abs(tk)); tk=tk(i); nav=nav(i,:);
A=nav(17)^2; e=nav(15); n=sqrt(MU/A^3)+nav(12); M=nav(13)+n*tk;
E=M; Ek=0; while abs(E-Ek)>1E-12, Ek=E; E=M+e*sin(Ek); end % solve kepler's eq.
p=atan2(sqrt(1-e^2)*sin(E),cos(E)-e)+nav(24);
r=[p;A*(1-e*cos(E));nav(22)+nav(26)*tk]+nav([16,14;11,23;21,19])*[sin(2*p);cos(2*p)];
OMG=nav(20)+(nav(25)-OMGE)*tk-OMGE*nav(18);
rs=r(2)*Rz(-OMG)*Rx(-r(3))*[cos(r(1));sin(r(1));0];
if nav(1)<70, y=2000; else y=1900; end
tc=(t0-datenum(nav(1)+y,nav(2),nav(3)))*86400+t-nav(4:6)*[3600;60;1];
dt=nav(7:9)*[1;tc;tc^2]-2*sqrt(MU*A)*e*sin(E)/C^2;