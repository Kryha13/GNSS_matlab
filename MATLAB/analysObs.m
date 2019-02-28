function [obsTable, obsMatrix, obsType] = analysObs(epoch, observablesHeader, xyzStation, obs, sysConst);
const;
noPRN = sysGNSS(sysConst);

signal = length(observablesHeader{sysConst})-5;
obsType = cellstr(observablesHeader{sysConst});
obsType = obsType(6:end);

%---obsTable: epoch, clocRec, signals: number of obs;
obsTable=repmat({nan}, length(epoch), signal+2);
obsTable(1:length(epoch), 1) = num2cell(epoch);

%---clockRec;
[~, ib, ic] = intersect(epoch, obs(:,2));
obsTable(ib, 2) = num2cell(obs(ic,3));

%---obsMatrix: observations (epoch x satellites x signals)
obsMatrix(1:length(epoch), 1:MaxNoPRN, 1:signal) = nan;


for i = 1:MaxNoPRN
        PRN = noPRN + i;
        obsSv = obs(obs(:,5) == PRN,:);
        [~, ie, ip] = intersect(epoch,obsSv(:,2));
        obsMatrix(ie, PRN-noPRN, :) = permute(obsSv(ip, 6:end),[1 3 2]);
end

for j = 1:length(epoch)
    for i = 1:signal
        %obsTable(j, i+2) = {find(~isnan(squeeze(obsMatrix(j,:,i))))};
        obsTable(j, i+2) = {find(~isnan(squeeze(obsMatrix(j,:,i))))+noPRN};
    end
end


%emptyEpoch = setdiff(epoch,obs(:,2)');
%squeeze 3wymiar na wektor
%cell2mat
%eval(['f', sysConst, type(2)]); string i char na zmiann¹