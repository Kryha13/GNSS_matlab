
function [XYZ_station,NEU_antenna,observablesHeader,measurementsInterval,timeFirst,timeLast,obsG,obsR,obsE,obsC,obsJ,obsI,obsS,slotR]=readRinex302(file)
%% This function open RINEX 3.02 observation files.
% Follows RINEX 3.02 standard. Reads Multiconstellation observables and
% generates an output matrix

%%% ------ Input--- %%%%%%%%%%
%
%   filePath : path to the RINEX observables file
%
%%% ------ Output--- %%%%%%%%%%
%
%   XYZ_station: ECEF coordinates of reference station or point (estimated by receiver)
%
%   observablesHeader: Cell array containing information of observables for
%   each constellation. to look for GPS information type
%   observablesHeader{'G'}
%
%   obs: Matrix containing observables {'week' 'epoch' 'flag' 'recClock' 'prn' 'C1C' 'D1C' 'S1C'}
%       Different for each constellation.
%
% This functions uses functions str2doubleq for quick conversion from
% string to double. It also uses function Date2GPSTime to convert a date to
% TOW.

const
%filePath = [folderObs,'\',file];

disp(['Reading rinex file : ',file]);
idFile  = fopen (file);

generalHeader = {'week', 'epoch', 'recClock', 'flag', 'prn'};

%Initialzie values
measurementsInterval = -1;
nObsSys(1:7)=nan;
obsG=[]; obsR=[]; obsE=[]; obsC=[]; obsJ=[]; obsI=[]; obsS=[]; 
slotR(1,1:MaxNoPRN)=nan;

%Progres indicator
nol = noLine(file);
nop = round(1:nol/10:nol);
nogl=0;
progress.text = sprintf( 'Completed: %%3.1f%%%%%%s' );

%% Read header
while (true)
    
    line = fgetl(idFile);                                                   %get line
    nogl=nogl+1;
    splitLine = strsplit(line);                                             %Line splited by spaces
    if isempty(splitLine{end})                                              %Delete last empty cell
        splitLine(end)=[];
    end        
    
    if (strfind(line,'APPROX POSITION XYZ'))                                % Receiver aprox position
       XYZ_station=[real(str2doubleq(splitLine(2))) real(str2doubleq(splitLine(3))) real(str2doubleq(splitLine(4)))];
    
    elseif (strfind(line,'ANTENNA: DELTA H/E/N'))                           % Antenna offsets
       NEU_antenna=[real(str2doubleq(splitLine(2))) real(str2doubleq(splitLine(3))) real(str2doubleq(splitLine(4)))]; 
    
    elseif (strfind(line,'TIME OF FIRST OBS'))
       timeFirst=[real(str2doubleq(splitLine(2))) real(str2doubleq(splitLine(3))) real(str2doubleq(splitLine(4))) real(str2doubleq(splitLine(5))) real(str2doubleq(splitLine(6))) real(str2doubleq(splitLine(7)))]; 

    elseif ~isempty(strfind(line,'SYS / # / OBS TYPES'))                    % Observation types for the different constellations (C1C, D1 and S1 only  )
        constellation = line(1);
        if constellation        == 'G'
            nObsSys(1) = real(str2doubleq(line(2:7)));
        elseif constellation    == 'R'
            nObsSys(2) = real(str2doubleq(line(2:7)));
        elseif constellation    == 'E'
            nObsSys(3) = real(str2doubleq(line(2:7)));
        elseif constellation    == 'C'
            nObsSys(4) = real(str2doubleq(line(2:7)));
        elseif constellation    == 'J'
            nObsSys(5) = real(str2doubleq(line(2:7)));
        elseif constellation    == 'I'
            nObsSys(6) = real(str2doubleq(line(2:7)));
        elseif constellation    == 'S'
            nObsSys(7) = real(str2doubleq(line(2:7)));
        end
        
      
        nObservables = str2doubleq(line(2:7));                              % Number of observables
        observables = splitLine(3:end - 6);                                 % Take the observables only (Not the text regions)
        observables = [generalHeader, observables];                         % Append the header, as the data will be ordered like specified in obsrvables now.
        
        if nObservables >13                                                 %Two line case
            line2 = fgetl(idFile);
            nogl=nogl+1;
            splitLine2 = strsplit(line2);
                if isempty(splitLine2{end})                                  %Delete last empty cell
                    splitLine2(end)=[];
                end
            observables = [observables, splitLine2(2:end - 6) ];
        end
        
        observablesHeader{uint8(constellation)} = observables;              % Contains all the observables for the constellations.
        %use constellation letter for indexation
        
    elseif strfind(line,'INTERVAL')
        measurementsInterval=real(str2doubleq(line(5:10)));                 % Measurement intervals (Default 1)
        
    elseif strfind(line,'GLONASS SLOT / FRQ #')
        nSlot = str2doubleq(line(1:4));                                     % Number of GLONASS satellites
        slot = splitLine(3:end - 5);                                        % Take the observables only (Not the text regions)
        if nSlot >8 %Two line case
            line2 = fgetl(idFile);
            nogl=nogl+1;
            splitLine2 = strsplit(line2);
            slot = [slot, splitLine2(2:end - 5) ];
        end
        for i = 1:real(nSlot)
            Rsv = slot{2*i-1};
            slotR(1,str2num(Rsv(2:end)))=str2num(slot{2*i});
        end
        
    elseif strfind(line,'END OF HEADER');
        break;                                                              % End of header loop
    end
end




if measurementsInterval == -1                                               %If itnerval not set interval = 1
    measurementsInterval = 1;
end
%% Read body

obs = [];                                                                   % Output matrix
epoch = 0;                                                                  % Epoch counter
nObs = 1;

while(~feof(idFile))                                                        % Until end of file
    line = fgetl(idFile);                                                   % Get line
    nogl=nogl+1;
   
    if ~isempty(find(nop==nogl, 1))
        fprintf( progress.text, nogl/nol*100, sprintf('\n') );
    end
    splitLine = strsplit(line);                                             % Split line by spaces
    
    if strfind(line, '>')                                                   % New epoch
        epoch = epoch + 1;
        %Read time
        year = str2doubleq(splitLine(2));
        month = str2doubleq(splitLine(3));
        day = str2doubleq(splitLine(4));
        hour = str2doubleq(splitLine(5));
        minute = str2doubleq(splitLine(6));
        second = str2doubleq(splitLine(7));
        time = [year, month, day, hour, minute, second];
        time=real(time);
        [gpsWeek,~,~,tow]=cal2gpstime(time(1),time(2),time(3),time(4),time(5),time(6)); %Transform date to seconds of week
                 
        currentEpoch = tow;                                                 % Output
        currentSatellites = str2doubleq(splitLine(9));                      % Satellite information
        currentFlag = str2doubleq(splitLine(8));                            % flag (use/not use)
        if length(splitLine)>9
            recClock = real(str2doubleq(splitLine(10)));                    % receiver clock correction (optional)
        else
            recClock = 0;
        end
        
    else
        error(['Loading not correct, satellites skiped : ',filePath])        % Simple check, it should never jump if using the right rinex version
                  
    end
    
    if currentSatellites == 0
        error(['No satellites in epoch : ',filePath])
        
    end
    
    for i = 1:real(currentSatellites)     % Read the epoch satellites
        line = fgetl(idFile);
        nogl=nogl+1;

        if ~isempty(find(nop==nogl, 1))
            fprintf( progress.text, nogl/nol*100, sprintf('\n') );
        end
 
        constellation = line(1);                                            % First character indicates the constellation
        prn = str2doubleq ([line(2) line(3)]);                              % Satellites PRN number
        
        nObservables = cellfun('length',observablesHeader(uint8(constellation))) - 3; %The header also includes things that are not measurements
        measurementsPosition = (4:16:16*nObservables+4);                    %Vector containing the columns of the measurements. Each 16 columns theres a measurement
        
        if measurementsPosition(end) > length(line)                         %Correction of a wierd bug
           measurementsPositionclear = find(measurementsPosition > length(line));
           measurementsPosition(measurementsPositionclear(2:end)) = []; 
           measurementsPosition(end) = length(line);       
           nObservables = length(measurementsPosition)-1;
        end
        
        measurementsValue = zeros(1,nObservables);                          %Initialize vector to store data

        for m = 1:nObservables                                              % Number of observables in the line (Generally 3)
            value = line(measurementsPosition(m):measurementsPosition(m+1));% Column position of measurement. Measurements take 16 columns
            if length(strsplit(value))>2                                    % Find the value if observation have quality number
               C =  strsplit(value);
               [~,b] = max(cellfun(@length, C));
               value = C(b);
            end
            measurementsValue(m) = str2doubleq(value);                      % convert string to double
        end
        
        measurementsValue = round(real(measurementsValue)*1000)/1000;       % output of str2doubleq is imaginary; round to 0.001
        if measurementsValue(1) == 0                                        % if PSR value equals 0
            continue;                                                       % Skip line (Satellite has no information on L1)
        end
        switch constellation                                                %Asign constellation based on first char of line
            
            case 'G' %GPS
                prn = prn+000;
            case 'R' %GLONASS
                prn = prn+100;
            case 'E' %Galileo
                prn = prn+200;
            case 'C' %BeiDou
                prn = prn+300;
            case 'J' %QZSS PRN-192
                prn = prn+1000;
            case 'I' %IRNSS
                prn = prn+2000;
            case 'S' %SBAS PRN-100
                prn = prn+3000;
            otherwise
                error(['Unrecognized constellation : ',filePath])
        end
        
        data = [gpsWeek,currentEpoch,recClock,currentFlag,prn,measurementsValue];   % store data
        obs{nObs} = real(data);
        nObs= nObs+1;
     
    end
  
  if ~isempty(find(nop==nogl, 1))
    fprintf( progress.text, nogl/nol*100, sprintf('\n') );
  end       
    
timeLast=time;

end
 
%Convert cell array to matrix. 

disp(['Convert to systems matrix : ',file]);
 
noG=0; noR=0; noE=0; noC=0; noJ=0; noI=0; noS=0;     
for i = 1:length(obs)
    Obs = cell2mat(obs(i));
    
    if Obs(5)>0 && Obs(5)<100
        noG=noG+1;
        if size(Obs,2) < nObsSys(1)+5;
            Obs(end+1:nObsSys(1)+5) = 0;
        end
        obsG(noG,:)=Obs;    
    
    elseif  Obs(5)>100 && Obs(5)<200
        noR=noR+1;
        if size(Obs,2) < nObsSys(2)+5;
            Obs(end+1:nObsSys(2)+5) = 0;
        end
        obsR(noR,:)=Obs;
    
    elseif  Obs(5)>200 && Obs(5)<300
        noE=noE+1;
        if size(Obs,2) < nObsSys(3)+5;
            Obs(end+1:nObsSys(3)+5) = 0;
        end
        obsE(noE,:)=Obs;
     
    elseif  Obs(5)>300 && Obs(5)<400
        noC=noC+1;
        if size(Obs,2) < nObsSys(4)+5;
            Obs(end+1:nObsSys(4)+5) = 0;
        end
        obsC(noC,:)=Obs;
    
    elseif  Obs(5)>1000 && Obs(5)<2000
        noJ=noJ+1;
        if size(Obs,2) < nObsSys(5)+5;
            Obs(end+1:nObsSys(5)+5) = 0;
        end
        obsJ(noJ,:)=Obs;
        
    elseif  Obs(5)>2000 && Obs(5)<3000
        noI=noI+1;
        if size(Obs,2) < nObsSys(6)+5;
            Obs(end+1:nObsSys(6)+5) = 0;
        end
        obsI(noI,:)=Obs;
        
    elseif  Obs(5)>3000 && Obs(5)<4000
        noS=noS+1;
        if size(Obs,2) < nObsSys(7)+5;
            Obs(end+1:nObsSys(7)+5) = 0;
        end
        obsS(noS,:)=Obs;
    end    
end

if ~isempty(obsG) temp=obsG(:,5:end); temp(temp==0)=nan; obsG=[obsG(:,1:4) temp]; end % obsG(obsG==0)=nan;
if ~isempty(obsR) temp=obsR(:,5:end); temp(temp==0)=nan; obsR=[obsR(:,1:4) temp]; end %obsR(obsR==0)=nan;
if ~isempty(obsE) temp=obsE(:,5:end); temp(temp==0)=nan; obsE=[obsE(:,1:4) temp]; end %obsE(obsE==0)=nan;
if ~isempty(obsC) temp=obsC(:,5:end); temp(temp==0)=nan; obsC=[obsC(:,1:4) temp]; end %obsC(obsC==0)=nan;
if ~isempty(obsJ) temp=obsJ(:,5:end); temp(temp==0)=nan; obsJ=[obsJ(:,1:4) temp]; end %obsJ(obsJ==0)=nan;
if ~isempty(obsI) temp=obsI(:,5:end); temp(temp==0)=nan; obsI=[obsI(:,1:4) temp]; end %obsI(obsI==0)=nan;
if ~isempty(obsS) temp=obsS(:,5:end); temp(temp==0)=nan; obsS=[obsS(:,1:4) temp]; end %obsS(obsS==0)=nan;

disp(['Observables loaded : ',file]);
fclose(idFile);


