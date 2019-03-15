
% file=['SEPT_20181740000_01H_01S_MO.rnx']; % 2018
% xyzStation = [3655333.1790  1403901.2258  5018038.3565];
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
% 

%% code position

% X = [3655333.1790  1403901.2258  5018038.3565; %SEPT
%     3655336.8160  1403899.0382  5018036.4527]; % SEP2
%  
% codes = ["C1C" "C5Q"]; % wybrane kody
time_interval = 100; % interwa³ na jaki chcemy pozycje
system = 'E';
% data = ["2018-SEPT.mat" "2018-SEP2.mat"]; %wczytane dane dla odbiorników
% 
% [X_code1, DOP1, dtrec1, Az_EL_R1, tropo1, tau1] = codepos(X(1,:), codes, time_interval, system, data(1));
% [X_code2, DOP2, dtrec2, Az_EL_R2, tropo2, tau2] = codepos(X(2,:), codes, time_interval, system, data(2));
% 
% d = X_code1 - X_code2

%% odl przybli¿ona po uwzglêdnieniu poprawki zegara odbiornika

% ro_approx1 = ro_approx(X_code1, dtrec1, tau1, data(1));
% ro_approx2 = ro_approx(X_code2, dtrec2, tau2, data(2));
% 
% dro = squeeze(ro_approx1(1,3,:)) - squeeze(Az_EL_R1(1,5,:)); %róznica miêdzy odl z kodowego a przybli¿on¹ dla odb1 sat1

%% pojedyncze ró¿nice

obs_types = ["C1C" "L1C" "C6C" "L6C" "C5Q" "L5Q" "C7Q" "L7Q" "C8Q" "L8Q"]; % wybór obserwacji do ró¿nic
data = 'data_all.mat';
[single_diff] = single_differences(time_interval, obs_types, system, obsMatrix1, obsMatrix2, Az_EL_R1, data);

%% podwójne róznice

double_diff = double_differences(single_diff, length(obs_types), 8);

