
% file=['SEP2_20181740000_01H_01S_MO.rnx']; % 2018
% xyzStation = [3655336.8160  1403899.0382  5018036.4527];
% measurementsInterval = 1;
% 
% [~, ~, observablesHeader, ~, timeFirst, timeLast, obsG, obsR, obsE, obsC, obsJ, obsI, obsS, slot]=readRinex302(file);
% [time(1), time(2), time(3), time(4)] = cal2gpstime(timeFirst);
% [~, ~, ~, gpsSecondsFirst] = cal2gpstime(timeFirst);
% [~, ~, ~, gpsSecondsLast] = cal2gpstime(timeLast);
% epoch = gpsSecondsFirst:measurementsInterval:gpsSecondsLast;
% sysDigit = [2]; % GPS=0, GLONASS=1, GALILEO=2, BEIDOU=3 --> choose systems
% 
% eph = readSp3(time, sysDigit);
% sysConst = ['E']; %inne systemy: 'G', 'R', 'C', 'E'
% sat_clk = readCLK(time, sysConst);
% bias = readBIA(time, sysConst);
% 
% 
% [obsTable, obsMatrix, obsType] = analysObs(epoch, observablesHeader, xyzStation, eval(['obs', sysConst]), sysConst);

X = [3655336.8160  1403899.0382  5018036.4527;
     3655333.1790  1403901.2258  5018038.3565]; % wsp z rinex
codes = ["C1C" "C5Q"]; % wybrane kody
time_interval = 1000; % interwa³ na jaki chcemy pozycje
system = 'E';

[X_code1, DOP1, dtrec1, Az_EL_R1] = codepos(X(1,:), codes, time_interval, system);
[X_code2, DOP2, dtrec2, Az_EL_R2] = codepos(X(2,:), codes, time_interval, system);

