
file=['SEPT_20173350000_01H_01S_MO.rnx'];
xyzStation = [3655333.908  1403901.088  5018038.137];
measurementsInterval = 1;

[~, ~, observablesHeader, ~, timeFirst, timeLast, obsG, obsR, obsE, obsC, obsJ, obsI, obsS, slot]=readRinex302(file);
[time(1), time(2), time(3), time(4)] = cal2gpstime(timeFirst);
[~, ~, ~, gpsSecondsFirst] = cal2gpstime(timeFirst);
[~, ~, ~, gpsSecondsLast] = cal2gpstime(timeLast);
epoch = gpsSecondsFirst:measurementsInterval:gpsSecondsLast;
eph = readSp3(time);

sysConst='E'; %inne systemy: 'G', 'R', 'C', 'E'
[obsTable, obsMatrix, obsType] = analysObs(epoch, observablesHeader, xyzStation, eval(['obs', sysConst]), sysConst);
     