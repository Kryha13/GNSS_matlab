% single point positioning -----------------------------------------------------
function [rr,t,dt,dop]=pointp(t0,obs,iobs,nav,inav,maska,xp)
C=299792458; f1=1.57542E9; f2=1.2276E9;

    i=find(iobs(:,3)==1); tr=iobs(i(1),1); sats=iobs(i,2);
    P(1:length(sats),1:length(sats))=0;
    %model troposfery------------------------------------------------
    [td,tw]=saasta(xp);
    for n=1:length(sats)
        [r,e,el(n),dts]=geodist(t0,tr,xp,0,nav(inav==sats(n),:));
        [md(n),mw(n)]=niell(t0,xp,el(n));
        trop(n)=td*md(n)+tw*mw(n);
        P(n,n)=1/cos(el(n));
    end
    obs(:,5)=obs(:,3)-trop';
    obs(:,6)=obs(:,4)-trop';
    %----------------------------------------------------------------
    y=obs(i,5:6)*[f1^2;-f2^2]/(f1^2-f2^2); % ion-free pseudorange
    %SPP-------------------------------------------------------------
    x=zeros(4,1); xk=ones(4,1);
    while norm(x-xk)>0.1  % dok³adnoœæ   
        [h,H]=prmodel(t0,tr,x(1:3),x(4)/C,sats,nav,inav,maska);       
        i=find(~isnan(y)&~isnan(h)); H=H(i,:);
        if length(i)<4, x(:)=nan; break, end
        xk=x; x=x+(H'*P(i,i)*H)\H'*P(i,i)*(y(i)-h(i));
    end
    rr=x(1:3); dt=x(4); t=tr-dt/C; % t=tr-dtr
    dop=sqrt(trace(inv(H'*H)));
