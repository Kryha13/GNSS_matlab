function [biases] = readBIA(fileTime, sysConst)

const
filePath = ['gbm', num2str(fileTime(1)), num2str(fileTime(2)), '.bia'];

disp(['Reading BIA file : ',filePath]);
idFile  = fopen (filePath);

bias = [];
nObs = 1;

while 1
    line=fgetl(idFile); if ~isstr(line), break, end 
    
    splitLine = strsplit(line);
            
            if (strfind(line(2:4), 'ISB'))
                    
                    if (strfind(line(7), sysConst))

                        constellation = cellstr(splitLine(3));
                        prn = splitLine(4);

                        station = splitLine(5);
                        obs1 = splitLine(6);
                        obs2 = splitLine(7);
                        bias_value = splitLine(11);
                        std = splitLine(12);

                        bias(nObs,:) = [constellation, prn, station, obs1, obs2, bias_value, std];
                        nObs= nObs+1;

                    end
            end
end
fclose('all');
end



