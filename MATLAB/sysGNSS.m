function [noPRN, nameSys] = sysGNSS(sysConst)

switch sysConst                                                
    case 'G' %GPS
        noPRN = 0;
        nameSys = 'GPS';
    case 'R' %GLONASS
        noPRN = 100;
        nameSys = 'GLONASS';
    case 'E' %Galileo
        noPRN = 200;
        nameSys = 'Galileo';
    case 'C' %BeiDou
        noPRN = 300;
        nameSys = 'BeiDou';
    case 'J' %QZSS PRN-192
        noPRN = 1000;
        nameSys = 'QZSS';
    case 'I' %IRNSS
        noPRN = 2000;
        nameSys = 'IRNSS';
    case 'S' %SBAS PRN-100
        noPRN = 3000;
        nameSys = 'SBAS';
end

end

