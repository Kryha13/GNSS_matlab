function[Ba,La,Ha]=hirvonen(X,Y,Z,a,e2);
j=1;
r = (X^2+Y^2)^(0.5);
B(j)=atan(Z/(r*(1-e2)));

while 1
   j = j+1; 
   N(j) = Np(B(j-1),a,e2);
   H(j) = (r/cos(B(j-1)))-N(j);
   B(j) = atan(Z/(r*(1-(e2*(N(j)/(N(j)+H(j)))))));
   
   if abs(B(j)-B(j-1))<(0.000001/206265) break; end;
end;

Ba = (B(j));
La = atan(Y/X);
Na = Np(Ba,a,e2);
Ha = (r/cos(Ba))-Na;

       