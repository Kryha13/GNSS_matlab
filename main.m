
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
% [obsTable, obsMatrix, obsType] = analysObs(epoch, observablesHeader, xyzStation, eval(['obs', sysConst]), sysConst);

%% pozycja kodowa

% X = [3655333.1790  1403901.2258  5018038.3565; %SEPT
%     3655336.8160  1403899.0382  5018036.4527]; % SEP2
%   
% codes = ["C1C" "C5Q"]; % wybrane kody
% time_interval = 30; % interwa³ na jaki chcemy pozycje
% system = 'E';
% data = ["2018-SEPT.mat" "2018-SEP2.mat"]; %wczytane dane dla odbiorników
% 
% [X_code1, DOP1, dtrec1] = codepos(X(1,:), codes, time_interval, system, data(1));
% [X_code2, DOP2, dtrec2] = codepos(X(2,:), codes, time_interval, system, data(2));

%%
% POWY¯SZE WCZYTANE I OBLICZONE ZAPISANE W PLIKU 'dane_do_test.mat'

%% Pozycja z obserwacji nieró¿nicowalnych

% time_interval = 30;
% obs_types = ["C1C" "L1C" "C6C" "L6C" "C7Q" "L7Q" "C5Q" "L5Q"];
% 
% [x_fixed, ratio] = calculate_position(obs_types, time_interval); 


%% Metoda CAR

[xN_step1, N_CAR] = CAR(obs_types, time_interval);




