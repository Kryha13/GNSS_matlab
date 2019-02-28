function [sat_clk] = readCLK(fileTime, sysConst)

const
filePath = ['gbm', num2str(fileTime(1)), num2str(fileTime(2)), '.clk'];

disp(['Reading CLK file : ',filePath]);
idFile  = fopen (filePath);

sat_clk = [];
nObs = 1;

while 1
    line=fgetl(idFile); if ~isstr(line), break, end 
    
    splitLine = strsplit(line);
    
%     for i=1:1000
        if (strfind(line(1:2),'AS'))
            
            for i=1:length(sysConst)
            
            if (strfind(line(4), sysConst(i)))
                
            constellation = line(4);
            prn = real(str2doubleq ([line(5) line(6)]));

            year = str2doubleq(splitLine(3));
            month = str2doubleq(splitLine(4));
            day = str2doubleq(splitLine(5));
            hour = str2doubleq(splitLine(6));
            minute = str2doubleq(splitLine(7));
            second = str2doubleq(splitLine(8));
            time = [year, month, day, hour, minute, second];
            time = real(time);
            [~, ~, ~, epoch] = cal2gpstime(time);

            sat_clk_corr = str2doubleq(splitLine(10));

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

%                 for i = 1:length(sysDigit)              %get data only for wanted system
%                     if round(prn/100)==sysDigit(i)
                        sat_clk(nObs,:) = [epoch, prn, sat_clk_corr];
                        nObs= nObs+1;
                  end
             end
        end
%     end
end
fclose('all');
end 

