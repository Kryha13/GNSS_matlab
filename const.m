%*******************************************************
%
% DESCRIPTION:
% 		This script contains many useful constants for GPS
%       and related work.  It should be kept in only
%       one place so that updates are immediately available 
%       to all other scripts/functions.
%  
%*******************************************************

% GENERAL CONSTANTS

% =========================================================================

% =========================================================================
MaxNoPRN = 40;
% =========================================================================
typeProd = 'gbm';        % sp3/clk/bia  -> gbm:GFZ, com:CODE
server = 'cddis.nasa.gov:21';
% =========================================================================
c = 299792458;           %----> Speed of light (meters/s).
Re = 6378137 ;           %----> Earth Radius (meters)
a = 6378137;             % GRS80 (meters)
e2 = 0.00669438002290;   % GRS80 
OMEGA = 7.2921151467E-5; % Earth rotation speed (rad/sec)
% =========================================================================
%%
fG1  = 1575420000;       %    L1-CARRIER FREQUENCY   GPS       1/SEC
fG2  = 1227600000;       %    L2-CARRIER FREQUENCY   GPS       1/SEC
fG5  = 1176450000;       %    L5-CARRIER FREQUENCY   GPS       1/SEC
%%
fR1  = 1602000000;       %    L1-CARRIER FREQUENCY   GLONASS   1/SEC
fR2  = 1246000000;       %    L2-CARRIER FREQUENCY   GLONASS   1/SEC
fR3  = 1202025000;       %    L3-CARRIER FREQUENCY   GLONASS   1/SEC
DfR1 =     562500;       %    L1-CARRIER FREQ. DIFF. GLONASS   1/SEC
DfR2 =     437500;       %    L2-CARRIER FREQ. DIFF. GLONASS   1/SEC
DfR3 =          0;       %    L3-CARRIER FREQ. DIFF. GLONASS   1/SEC
%%
fE1  = 1575420000;       %    L1-CARRIER FREQUENCY   GALILEO   1/SEC
fE8  = 1191795000;       %    L5-CARRIER FREQUENCY   GALILEO   1/SEC
fE5  = 1176450000;       %    L5a-CARRIER FREQUENCY  GALILEO   1/SEC
fE7  = 1207140000;       %    L5b-CARRIER FREQUENCY  GALILEO   1/SEC
fE6  = 1278750000;       %    L6-CARRIER FREQUENCY   GALILEO   1/SEC
%%
fS1  = 1575420000;       %    L1-CARRIER FREQUENCY   SBAS      1/SEC
fS5  = 1176450000;       %    L5-CARRIER FREQUENCY   SBAS      1/SEC
fC1  = 1589740000;       %    L1-CARRIER FREQUENCY   COMPASS   1/SEC
fC2  = 1561098000;       %    L2-CARRIER FREQUENCY   COMPASS   1/SEC
fC7  = 1207140000;       %    L5b-CARRIER FREQUENCY  COMPASS   1/SEC
fC6  = 1268520000;       %    L6-CARRIER FREQUENCY   COMPASS   1/SEC
fJ1  = 1575420000;       %    L1-CARRIER FREQUENCY   QZSS      1/SEC
fJ2  = 1227600000;       %    L2-CARRIER FREQUENCY   QZSS      1/SEC
fJ5  = 1176450000;       %    L5-CARRIER FREQUENCY   QZSS      1/SEC
fJ6  = 1278750000;       %    L6-CARRIER FREQUENCY   QZSS      1/SEC
fI5  = 1176450000;       %    L5-CARRIER FREQUENCY   IRNSS     1/SEC
fI9  = 2492028000;       %    S-CARRIER FREQUENCY    IRNSS     1/SEC
% =========================================================================
obsT = {'C', 'L', 'S'};
obsTypeG = {'1C','1W','2L','2W','5Q'};
obsTypeR = {'1C','1P','2C','2P','3Q'};
obsTypeE = {'1C','5Q','7Q','8Q','6C'};
obsTypeC = {'2I','7I','6I'};
obsTypeI = {'5A'};
obsTypeS = {'1C','5I'};
obsTypeJ = {'1C','1X','1Z','2X','5X'};
% =========================================================================
noiseOutliersC = 6;           % (meters)
noiseOutliersL = 0.05;        % (meters) 
noiseOutliersSNR = 5;        % (dB-Hz) 
noiseOutliersMethod = 'mean';       % outliers detection method if MATLAB (function isoutlier)
                         % others: 'median','mean','quartiles','grubbs','gesd'

%{
GM     = 398.6004415D12   GRAVITY CONSTANT*EARTH MASS      M**3/SEC**2
GMS    = 1.3271250D20     GRAVITY CONSTANT*SOLAR MASS      M**3/SEC**2
GMM    = 4.9027890D12     GRAVITY CONSTANT*LUNAR MASS      M**3/SEC**2
AU     = 149597870691     ASTRONOMICAL UNIT                M
AE     =    6378137.D0    EQUATORIAL RADIUS OF EARTH       M
CONRE  =    6371000.D0    MEAN RADIUS OF THE EARTH         M
J2     =   1.0826359D-3   DYNAMICAL FORM-FACTOR IERS(2003) 1
FACTEC = 40.3D16          IONOSPHERIC FACTOR               M/SEC**2/TECU
P0     =          -.94D-7 NOMINAL RAD.PR. ACCELERAT.       M/SEC**2
OMEGA  = 7292115.1467D-11 ANGULAR VELOCITY OF EARTH        RAD/SEC
ET-UTC =         55.      EPH. TIME (ET) MINUS UTC         SEC
WGTPHA =          1.D0    WEIGHT FOR PHASE OBSERVATIONS    1
WGTCOD =          1.D-4   WEIGHT FOR CODE OBSERVATIONS     1
HREF   =          0.      REFERENCE HEIGHT FOR METEO MODEL M
PREF   =       1013.25    PRESSURE AT HREF                 MBAR
TREF   =         18.      TEMPERATURE AT HREF              DEG. CELSIUS
HUMREF =         50.      HUMIDITY AT HREF                 %
%}
% =========================================================================
