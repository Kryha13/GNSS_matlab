function [pos, dtsat]=wspsat(sv, eph, czas)
icol = find_eph(eph, sv, czas);
pos  = satpos (czas, eph(:, icol));
toc = eph(21,icol);
dt = check_t(czas-toc);
dtsat = (eph(2,icol)*dt + eph(20,icol))*dt + eph(19,icol);

function satp = satpos(t,eph);
%SATPOS   Calculation of X,Y,Z coordinates at time t
%                   for given ephemeris eph

%Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%$Revision: 1.1 $  $Date: 2004/02/09  $

GM = 3.986005e14;             % earth's universal gravitational
% parameter m^3/s^2
Omegae_dot = 7.2921151467e-5; % earth rotation rate, rad/s

%  Units are either seconds, meters, or radians
%  Assigning the local variables to eph
svprn   =   eph(1);
af2     =   eph(2);
M0      =   eph(3);
roota   =   eph(4);
deltan  =   eph(5);
ecc     =   eph(6);
omega   =   eph(7);
cuc     =   eph(8);
cus     =   eph(9);
crc     =  eph(10);
crs     =  eph(11);
i0      =  eph(12);
idot    =  eph(13);
cic     =  eph(14);
cis     =  eph(15);
Omega0  =  eph(16);
Omegadot=  eph(17);
toe     =  eph(18);
af0     =  eph(19);
af1     =  eph(20);
toc     =  eph(21);

% Procedure for coordinate calculation
A = roota*roota;
tk = check_t(t-toe);
n0 = sqrt(GM/A^3);
n = n0+deltan;
M = M0+n*tk;
M = rem(M+2*pi,2*pi);
E = M;
for i = 1:10
   E_old = E;
   E = M+ecc*sin(E);
   dE = rem(E-E_old,2*pi);
   if abs(dE) < 1.e-12
      break;
   end
end
E = rem(E+2*pi,2*pi);
v = atan2(sqrt(1-ecc^2)*sin(E), cos(E)-ecc);
phi = v+omega;
phi = rem(phi,2*pi);
u = phi              + cuc*cos(2*phi)+cus*sin(2*phi);
r = A*(1-ecc*cos(E)) + crc*cos(2*phi)+crs*sin(2*phi);
i = i0+idot*tk       + cic*cos(2*phi)+cis*sin(2*phi);
Omega = Omega0+(Omegadot-Omegae_dot)*tk-Omegae_dot*toe;
Omega = rem(Omega+2*pi,2*pi);
x1 = cos(u)*r;
y1 = sin(u)*r;
satp(1,1) = x1*cos(Omega)-y1*cos(i)*sin(Omega);
satp(2,1) = x1*sin(Omega)+y1*cos(i)*cos(Omega);
satp(3,1) = y1*sin(i);
%%%%%%%%% end satpos.m %%%%%%%%%

function tt = check_t(t);
% CHECK_T repairs over- and underflow of GPS time

%  Written by Kai Borre
%  April 1, 1996

 half_week = 302400;
 tt = t;

 if t >  half_week, tt = t-2*half_week; end
 if t < -half_week, tt = t+2*half_week; end

%%%%%%% end check_t.m  %%%%%%%%%%%%%%%%%

function icol = find_eph(Eph,sv,time)
%FIND_EPH  Finds the proper column in ephemeris array

%Kai Borre and C.C. Goad 11-26-96
%Copyright (c) by Kai Borre
%$Revision: 1.2 $  $Date: 2004/02/09  $

icol = 0;
isat = find(Eph(1,:) == sv);
n = size(isat,2);
if n == 0
   return
end;
icol = isat(1);
dtmin = Eph(21,icol)-time;
for t = isat
   dt = Eph(21,t)-time;
   if dt < 0
      if abs(dt) < abs(dtmin)
         icol = t;
         dtmin = dt;
      end
   end
end
%%%%%%%%%%%%  find_eph.m  %%%%%%%%%%%%%%%%%
