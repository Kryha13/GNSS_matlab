% pseudorange model (ionosphere-free) ------------------------------------------
function [h,H]=prmodel(t0,tr,rr,dtr,sats,nav,inav,maska)
%t0: data
%tr: epoka obserwacji
%rr: wektor pozycji
%dtr: dt/c [s]
%sats: PRN SV dla danej epoki

C=299792458;
for n=1:length(sats)
    [r,e,el(n),dts]=geodist(t0,tr-dtr,rr,maska,nav(inav==sats(n),:));
    % r:odleg³oœæ geometryczna obd-sat (d³ugoœæ)
    % e:wektor jednostkowy rrs
    % el(n):elewacja satelity
    % dts:poprawka zegara sat. 
   
    h(n,1)=r+C*(dtr-dts); H(n,:)=[e',1];
end