function v=vpr(t0,tr,dt,obs,iobs,nav,inav,x)
C=299792458; f1=1.57542E9; f2=1.2276E9;
i=find(iobs(:,3)==1); sats=iobs(i,2);

[td,tw]=saasta(x);
for n=1:length(sats)
    [r(n,1),e,el(n),dts(n)]=geodist(t0,tr-dt/C,x,0,nav(inav==sats(n),:));
    [md(n),mw(n)]=niell(t0,x,el(n));
    trop(n)=td*md(n)+tw*mw(n);
end

obs(:,5)=obs(:,3)-trop'+dts'*C-dt;
obs(:,6)=obs(:,4)-trop'+dts'*C-dt;
i=find(~isnan(obs(:,5))&~isnan(r));
v(:,1)=r(i)-obs(i,5);
v(:,2)=r(i)-obs(i,6);
v(:,3)=sats(i);

