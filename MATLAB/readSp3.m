function [obs] = readSp3(fileTime)

const
filePath = ['gbm', num2str(fileTime(1)), num2str(fileTime(2)), '.sp3'];

disp(['Reading SP3 file : ',filePath]);
idFile  = fopen (filePath);

obs = [];
nObs = 1;

while 1
    line=fgetl(idFile); if ~isstr(line), break, end
    splitLine = strsplit(line);  
    
    if (strfind(line,'##')) line=fgetl(idFile); nSV = real(str2doubleq(line(4:7))); end
    
    if (strfind(line(1:1),'*'))
        year = str2doubleq(splitLine(2));
        month = str2doubleq(splitLine(3));
        day = str2doubleq(splitLine(4));
        hour = str2doubleq(splitLine(5));
        minute = str2doubleq(splitLine(6));
        second = str2doubleq(splitLine(7));
        time = [year, month, day, hour, minute, second];
        time = real(time);
        [~, ~, ~, epoch] = cal2gpstime(time);
        
        for i = 1:nSV    
            line = fgetl(idFile);
            splitLine = strsplit(line); 
            
            constellation = line(2);                                            % First character indicates de constellation
            prn = real(str2doubleq ([line(3) line(4)])); 
            X = real(str2doubleq(splitLine(2)))*1000; %[m]
            Y = real(str2doubleq(splitLine(3)))*1000;
            Z = real(str2doubleq(splitLine(4)))*1000;
            dt = real(str2doubleq(splitLine(5)))/1000000; %[s]
            
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
        
            obs(nObs,:) = [epoch,prn,X,Y,Z,dt];
            nObs= nObs+1;
        end            
    end
end
fclose('all');

end

