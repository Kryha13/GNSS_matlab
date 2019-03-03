function [obs,iobs,nav,inav]=readrnx(t0,ts,te,files)
% READRNX : read rinex observation data/navigation message files
%
% Read RINEX (Receiver Independent Exchange Format) GPS/GNSS observation data
% (OBS) and navigation message (NAV) files. Only supports GPS.
%
%          Copyright (C) 2006 by T.TAKASU, All rights reserved.
%
% argin  : t0    = date number by datenum.m (days)
%          ts,te = sampling start/end time relative to t0 (sec)
%          files = RINEX OBS/NAV file paths (cell array)
%
% argout : obs   = observation data sorted by sampling time (nan:no data)
%              obs(n,1)  : L1 phase (cycle)
%              obs(n,2)  : L2 phase (cycle)
%              obs(n,3)  : C1 pseudorange (m)
%              obs(n,4)  : P2 pseudorange (m)
%          iobs  = observation data index
%              iobs(n,1) : obs(n,:) sampling time relative to t0 (sec)
%              iobs(n,2) : obs(n,:) satellite PRN number
%              iobs(n,3) : obs(n,:) receiver index (count up file by file)
%          nav   = navigation messages
%              nav(n,1:6): Toc [year(2digits),month,day,hour,min,sec]
%              nav(n,7:9): SV Clock [bias,drift,drift-rate]
%              nav(n,10) : IODE,    nav(n,11) : Crs,     nav(n,12) : Delta n
%              nav(n,13) : M0,      nav(n,14) : Cuc,     nav(n,15) : e
%              nav(n,16) : Cus,     nav(n,17) : sqrt(A), nav(n,18) : Toe
%              nav(n,19) : Cic,     nav(n,20) : OMEGA0,  nav(n,21) : Cis
%              nav(n,22) : i0,      nav(n,23) : Crc,     nav(n,24) : omega
%              nav(n,25) : OMEGADOT,nav(n,26) : IDOT,    nav(n,27) : Codes on L2
%              nav(n,28) : GPS Week #,    nav(n,29) : L2 P data flag
%              nav(n,30) : SV accuracy,   nav(s,31) : SV health
%              nav(s,32) : TGD,           nav(n,33) : IODC
%              nav(n,34) : Trans. time,   nav(n,35:37): spare
%          inav  = navigation messages index
%              inav(n,1) : nav(n,:) satellite PRN number
%
% version: $Revision: 2 $ $Date: 06/01/29 7:15 $
% history: 2006/01/26 1.1 new

obs=[]; iobs=[]; nav=[]; inav=[]; rcv=0;
for n=1:length(files)
    f=fopen(files{n},'rt');
    if f<0, error(['file open error :',files{n}]); end
    disp(['reading rinex file ... : ',files{n}]);
    [trnx,tobs]=readrnxh(f);
    if trnx=='O'
        [o,i]=readrnxo(t0,ts,te,f,tobs);
        rcv=rcv+1; i(:,3)=rcv; obs=[obs;o]; iobs=[iobs;i];
    elseif trnx=='N'
        [v,i]=readrnxn(f); nav=[nav;v]; inav=[inav;i];
    end
    fclose(f);
end
if ~isempty(obs), [iobs,i]=sortrows(iobs); obs=obs(i,:); end
if ~isempty(nav)
    [n,i]=unique(nav(:,1:33),'rows'); nav=nav(i,:); inav=inav(i,:);
end

% read rinex header -----------------------------------------------------------
function [trnx,tobs]=readrnxh(f)
trnx=''; tobs={};
while 1
    s=fgetl(f); if ~isstr(s), break, end
    label=subs(s,61,20);
    if findstr(label,'RINEX VERSION / TYPE')
        trnx=subs(s,21,1);
    elseif findstr(label,'# / TYPES OF OBSERV')
        p=11;
        for i=1:s2n(s,1,6)
            if p>=59, s=fgetl(f); p=11; end
            tobs{i}=subs(s,p,2); p=p+6;
        end
    elseif findstr(label,'END OF HEADER'), break, end
end

% read rinex observation data --------------------------------------------------
function [obs,iobs]=readrnxo(t0,ts,te,f,tobs)
to={'L1','L2','C1','P2'}; ind=zeros(1,length(tobs));
for n=1:4, i=find(strcmp(to{n},tobs)); if ~isempty(i), ind(i)=n; end, end
obs=repmat(nan,100000,4); iobs=zeros(100000,3); n=0;
while 1
    s=fgetl(f); if ~isstr(s), break, end
    e=s2e(s,1,26);
    if length(e)>=6
        if e(1)<80, y=2000; else y=1900; end
        tt=(datenum(e(1)+y,e(2),e(3))-t0)*86400+e(4:6)*[3600;60;1];
        if round(tt)>te, break, else flag=round(tt)>=ts; end
        ns=s2n(s,30,3); iobs(n+(1:ns),1)=tt; p=33;
        for i=1:ns
            if p>=69, s=fgetl(f); p=33; end
            if flag&any(subs(s,p,1)==' G'), iobs(n+i,2)=s2n(s,p+1,2); end, p=p+3;
        end
        for i=1:ns
            s=fgetl(f); p=1;
            for j=1:length(tobs)
                if p>=80, s=fgetl(f); p=1; end
                if flag&ind(j)>0, obs(n+i,ind(j))=s2n(s,p,14); end, p=p+16;
            end
        end
        if flag, n=n+ns; end
    end
end
obs=obs(1:n,:); obs(obs==0)=nan; iobs=iobs(1:n,:); 

% read rinex navigation message ------------------------------------------------
function [nav,inav]=readrnxn(f)
nav=zeros(2000,37); inav=zeros(2000,1); n=0; m=1;
while 1
    s=fgetl(f); if ~isstr(s) break, end
    s=strrep(s,'D','E');
    if m==1, prn=s2n(s,1,2); a=s2e(s,3,30); else a=[a,s2n(s,4,19)]; end
    for p=23:19:61, a=[a,s2n(s,p,19)]; end
    if m<8, m=m+1; else n=n+1; nav(n,:)=a; inav(n,1)=prn; m=1; end
end
nav=nav(1:n,:); inav=inav(1:n,:); 

% string to number/epoch/substring  --------------------------------------------
function a=s2n(s,p,n), a=sscanf(subs(s,p,n),'%f'); if isempty(a), a=nan; end
function e=s2e(s,p,n), e=sscanf(subs(s,p,n),'%d %d %d %d %d %f',6)';
function s=subs(s,p,n), s=s(p:min(p+n-1,length(s)));
