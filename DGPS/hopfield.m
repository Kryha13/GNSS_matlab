function [td,tw]=hopfield(p,temp,hum);

c1=77.624;   c2=-12.924;  c3=3.719*10^5;
Nd0=c1*p/temp;
Nw0=(c2*hum/temp)+(c3*hum/temp^2);
hd=40136+148.72*(temp-273.15);
hw=11000;
%md=1/sin(deg2rad((el^2+6.25)^(0.5)));
%mw=1/sin(deg2rad((el^2+2.25)^(0.5)));
%trop=(10^-6/5)*(Nd0*hd*md+Nw0*hw*mw);
td=(10^-6/5)*Nd0*hd;
tw=(10^-6/5)*Nw0*hw;