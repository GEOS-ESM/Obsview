% GSIINFO: GSI-specific information.

OBSATTRIBUTES.names = {'kx','ks','kt','time','lat','lon','lev',...
    'obs','sigo','omf','oma','qch','qcx','xm'};
OBSATTRIBUTES.descr = {'data source','sounding index','data type',...
    'observation time','latitude','longitude','level/channel',...
    'observed value','observation error',...
    'background residual','analysis residual',...
    'QC history','QC exclusion','metadata'};
OBSATTRIBUTES.discr = [1 1 1 0 0 0 0 0 0 0 0 1 1 0];

% datatype descriptions and units:
% NOTE: data with sigo>msigo are considered excluded; see odsclean
% -------------------------------
KTS(1).value = 44;
KTS(1).id = 'Upper-air virtual temperature';
KTS(1).units = 'Kelvin';
KTS(1).msigo = 10;
KTS(2).value = 33;
KTS(2).id = 'Surface (2m) pressure';
KTS(2).units = 'hPa';
KTS(2).msigo = 12;
KTS(3).value = 4;
KTS(3).id = 'Upper-air zonal wind';
KTS(3).units = 'm/sec';
KTS(3).msigo = 20;
KTS(4).value = 5;
KTS(4).id = 'Upper-air meridional wind';
KTS(4).units = 'm/sec';
KTS(4).msigo = 20;
KTS(5).value = 11;
KTS(5).id = 'Upper-air specific humidity';
KTS(5).units = 'g/kg';
KTS(5).msigo = 20;
KTS(6).value = 12;
KTS(6).id = 'Surface (10m) wind speed';
KTS(6).units = 'm/sec';
KTS(6).msigo = 20;
KTS(7).value = 40;
KTS(7).id = 'Brightness temperature';
KTS(7).units = 'Kelvin';
KTS(7).msigo = 20;
KTS(8).value = 39;
KTS(8).id = 'Sea-surface temperature';
KTS(8).units = 'Kelvin';
KTS(8).msigo = 10;
KTS(9).value = 21;
KTS(9).id = 'Total column ozone';
KTS(9).units = 'Dobson Units';
KTS(9).msigo = 10;
KTS(10).value = 22;
KTS(10).id = 'Layer ozone';
KTS(10).units = 'Dobson Units';
KTS(10).msigo = 10;

% pressure-level data ('lev' attribute = pressure):
% -------------------
KTPRS = [4 5 11 44];

% surface data ('lev' attribute not used):
% ------------
KTSFC = [12 33 39];

% radiance data ('lev' attribute = channel number):
% -------------
KTRAD = [40];

% wind vectors
% ------------
KTWND = {[4 5]};

% data source descriptions
% -------------------
% (from http://www.emc.ncep.noaa.gov/mmb/papers/keyser/prepbufr.doc/table_2.htm)
KXS( 1).value = 102;
KXS( 1).id    = 'SSM/I 7-CHANNEL BRIGHTNESS TEMPERATURES (DMSP-13, DMSP-15)';
KXS( 2).value = 111;
KXS( 2).id    = 'SYNTHETIC TROPICAL CYCLONE STORM CENTER PSFC, Q';
KXS( 3).value = 120;
KXS( 3).id    = 'RAWINSONDE VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE';
KXS( 4).value = 122;
KXS( 4).id    = 'CLASS SOUNDING VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE';
KXS( 5).value = 130;
KXS( 5).id    = 'AIREP AND PIREP AIRCRAFT SENSIBLE TEMPERATURE';
KXS( 6).value = 131;
KXS( 6).id    = 'ASDAR AIRCRAFT SENSIBLE TEMPERATURE';
KXS( 7).value = 132;
KXS( 7).id    = 'FLIGHT-LEVEL RECONNAISSANCE AND PROFILE DROPSONDE VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE';
KXS( 8).value = 133;
KXS( 8).id    = 'MDCARS AIRCRAFT SENSIBLE TEMPERATURE (SPECIFIC HUMIDITY FLAGGED FOR NON-USE BY ANALYSIS)';
KXS( 9).value = 141;
KXS( 9).id    = 'INDIA IR AND VISIBLE CLOUD TEMPERATURE (INSAT-2E)';
KXS(10).value = 142;
KXS(10).id    = 'JMA IR AND VISIBLE CLOUD TEMPERATURE (GMS-5)';
KXS(11).value = 143;
KXS(11).id    = 'EUMETSAT IR AND VISIBLE CLOUD TEMPERATURE (METEOSAT-5, METEOSAT-7)';
KXS(12).value = 144;
KXS(12).id    = 'NESDIS VISIBLE CLOUD TEMPERATURE (GOES-8, GOES-10)';
KXS(13).value = 145;
KXS(13).id    = 'NESDIS IR CLOUD TEMPERATURE (GOES-8, GOES-10)';
KXS(14).value = 146;
KXS(14).id    = 'NESDIS IMAGER WATER VAPOR CLOUD TEMPERATURE AT CLOUD TOP (GOES-8, GOES-10)';
KXS(15).value = 147;
KXS(15).id    = 'NESDIS IMAGER WATER VAPOR CLOUD TEMPERATURE - DEEP LAYER (GOES-8, GOES-10)';
KXS(16).value = 148;
KXS(16).id    = 'NESDIS SOUNDER WATER VAPOR CLOUD TEMPERATURE AT CLOUD TOP (GOES-8, GOES-10)';
KXS(17).value = 149;
KXS(17).id    = 'NESDIS SOUNDER WATER VAPOR CLOUD TEMPERATURE - DEEP LAYER (GOES-8, GOES-10)';
KXS(18).value = 150;
KXS(18).id    = 'SSM/I SUPEROBED FNOC RAIN RATE (DMSP-13, DMSP-15)';
KXS(19).value = 151;
KXS(19).id    = 'NESDIS SFOV CLOUD TOP PRESSURE AND TEMPERATURE, CLOUD AMOUNT (GOES-8, GOES-10)';
KXS(20).value = 152;
KXS(20).id    = 'SSM/I SUPEROBED NEURAL NET 3 TOTAL PRECIPITABLE WATER RETRIEVALS (DMSP-13, DMSP-15)';
KXS(21).value = 156;
KXS(21).id    = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER LAND - CLEAR (GOES-8, GOES-10)';
KXS(22).value = 157;
KXS(22).id    = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER LAND - CLOUD CORRECTED (GOES-8, GOES-10)';
KXS(23).value = 158;
KXS(23).id    = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER OCEAN - CLEAR (GOES-8, GOES-10)';
KXS(24).value = 159;
KXS(24).id    = 'NESDIS LAYER PRECIPITABLE WATER RETRIEVALS OVER OCEAN - CLOUD CORRECTED (GOES-8, GOES-10)';
KXS(25).value = 164;
KXS(25).id    = 'NESDIS RADIANCES OVER LAND - CLEAR (GOES-8, GOES-10)';
KXS(26).value = 165;
KXS(26).id    = 'NESDIS RADIANCES OVER LAND - CLOUD CORRECTED (GOES-8, GOES-10)';
KXS(27).value = 174;
KXS(27).id    = 'NESDIS RADIANCES OVER OCEAN - CLEAR (GOES-8, GOES-10)';
KXS(28).value = 175;
KXS(28).id    = 'NESDIS RADIANCES OVER OCEAN - CLOUD CORRECTED (GOES-8, GOES-10)';
KXS(29).value = 180;
KXS(29).id    = 'SURFACE MARINE (SHIP, BUOY, C-MAN) VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE (STATION PRESSURE REPORTED)';
KXS(30).value = 181;
KXS(30).id    = 'SURFACE LAND SYNOPTIC AND METAR STATION PRESSURE, SPECIFIC HUMIDITY (TEMPERATURE NOT USED BY ANALYSIS) (STATION PRESSURE REPORTED)';
KXS(31).value = 182;
KXS(31).id    = 'SPLASH LEVEL VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE (OVER OCEAN ONLY)';
KXS(32).value = 183;
KXS(32).id    = 'SURFACE MARINE (SHIP, BUOY, C-MAN), LAND SYNOPTIC AND METAR VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE (STATION PRESSURE NOT REPORTED)';
KXS(33).value = 187;
KXS(33).id    = 'SURFACE METAR VIRTUAL TEMPERATURE, SPECIFIC HUMIDITY, STATION PRESSURE (STATION PRESSURE NOT REPORTED)';
KXS(34).value = 190;
KXS(34).id    = 'OPC/NOS POINT MEAN SEA-LEVEL PRESSURE BOGUS';
KXS(35).value = 191;
KXS(35).id    = 'AUSTRALIAN PAOB MEAN SEA-LEVEL PRESSURE BOGUS';
KXS(36).value = 210;
KXS(36).id    = 'SYNTHETIC TROPICAL CYCLONE U, V';
KXS(37).value = 220;
KXS(37).id    = 'RAWINSONDE U, V';
KXS(38).value = 221;
KXS(38).id    = 'PIBAL U, V';
KXS(39).value = 222;
KXS(39).id    = 'CLASS SOUNDING U, V';
KXS(40).value = 223;
KXS(40).id    = 'WIND PROFILER U, V';
KXS(41).value = 224;
KXS(41).id    = 'NEXRAD VERTICAL AZIMUTH DISPLAY (VAD) U, V';
KXS(42).value = 225;
KXS(42).id    = 'NEXRAD RADIAL U, V';
KXS(43).value = 230;
KXS(43).id    = 'AIREP AND PIREP AIRCRAFT U, V';
KXS(44).value = 231;
KXS(44).id    = 'ASDAR AIRCRAFT U, V';
KXS(45).value = 232;
KXS(45).id    = 'FLIGHT-LEVEL RECONNAISSANCE AND PROFILE DROPSONDE U, V';
KXS(46).value = 233;
KXS(46).id    = 'MDCARS AIRCRAFT U, V';
KXS(47).value = 241;
KXS(47).id    = 'INDIA IR AND VISIBLE CLOUD DRIFT U, V (INSAT-2E)';
KXS(48).value = 242;
KXS(48).id    = 'JMA IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS BELOW 850 MB (GMS-5)';
KXS(49).value = 243;
KXS(49).id    = 'EUMETSAT IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS BELOW 850 MB (METEOSAT-5, METEOSAT-7)';
KXS(50).value = 245;
KXS(50).id    = 'NESDIS IR CLOUD DRIFT U, V (GOES-8, GOES-10)';
KXS(51).value = 246;
KXS(51).id    = 'NESDIS IMAGER WATER VAPOR CLOUD U, V AT CLOUD TOP (GOES-8, GOES-10)';
KXS(52).value = 247;
KXS(52).id    = 'NESDIS IMAGER WATER VAPOR CLOUD U, V - DEEP LAYER (GOES-8, GOES-10)';
KXS(53).value = 248;
KXS(53).id    = 'NESDIS SOUNDER WATER VAPOR CLOUD U, V AT CLOUD TOP (GOES-8, GOES-10)';
KXS(54).value = 249;
KXS(54).id    = 'NESDIS SOUNDER WATER VAPOR CLOUD U, V - DEEP LAYER (GOES-8, GOES-10)';
KXS(55).value = 250;
KXS(55).id    = 'JMA WATER VAPOR CLOUD U, V (GMS-5)';
KXS(56).value = 251;
KXS(56).id    = 'NESDIS VISIBLE CLOUD DRIFT U, V (GOES-8, GOES-10)';
KXS(57).value = 252;
KXS(57).id    = 'JMA IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS ABOVE 850 MB (GMS-5)';
KXS(58).value = 253;
KXS(58).id    = 'EUMETSAT IR AND VISIBLE CLOUD DRIFT U, V AT LEVELS ABOVE 850 MB (METEOSAT-5, METEOSAT-7)';
KXS(59).value = 254;
KXS(59).id    = 'EUMETSAT WATER VAPOR CLOUD U, V (METEOSAT-5, METEOSAT-7)';
KXS(60).value = 255;
KXS(60).id    = 'NESDIS PICTURE TRIPLET CLOUD U, V (GOES-8, GOES-10)';
KXS(61).value = 256;
KXS(61).id    = 'INDIA WATER VAPOR CLOUD U, V (INSAT-2E)';
KXS(62).value = 280;
KXS(62).id    = 'SURFACE MARINE (SHIP, BUOY, C-MAN) U, V (STATION PRESSURE REPORTED)';
KXS(63).value = 281;
KXS(63).id    = 'SURFACE LAND SYNOPTIC AND METAR U, V (STATION PRESSURE REPORTED)';
KXS(64).value = 282;
KXS(64).id    = 'ATLAS BUOY U, V';
KXS(65).value = 283;
KXS(65).id    = 'SSM/I SUPEROBED NEURAL NET 3 WIND SPEED (DIRECTION SET TO ZERO, SPEED ASSIMILATED DIRECTLY) (DMSP-13, DMSP-15)';
KXS(66).value = 284;
KXS(66).id    = 'SURFACE MARINE (SHIP, BUOY, C-MAN), LAND SYNOPTIC AND METAR U, V (STATION PRESSURE NOT REPORTED)';
KXS(67).value = 285;
KXS(67).id    = 'QUIKSCAT SCATTEROMETER U, V';
KXS(68).value = 286;
KXS(68).id    = 'ERS-2 SCATTEROMETER U, V';
KXS(69).value = 287;
KXS(69).id    = 'SURFACE METAR U, V (STATION PRESSURE NOT REPORTED)';

KXS(70).value = 14;
KXS(70).id    = 'NOAA-14 HIRS/2';
KXS(71).value = 214;
KXS(71).id    = 'NOAA-14 MSU';
KXS(72).value = 315;
KXS(72).id    = 'NOAA-15 AMSUA';
KXS(73).value = 415;
KXS(73).id    = 'NOAA-15 AMSUB';
KXS(74).value = 60;
KXS(74).id    = 'GOES-10 Sounder';
KXS(75).value = 62;
KXS(75).id    = 'GOES-12 Sounder';
KXS(76).value = 316;
KXS(76).id    = 'NOAA-16 AMSUA';
KXS(77).value = 416;
KXS(77).id    = 'NOAA-16 AMSUB';
KXS(78).value = 16;
KXS(78).id    = 'NOAA-16 HIRS/3';
KXS(79).value = 417;
KXS(79).id    = 'NOAA-17 AMSUB';
KXS(80).value = 17;
KXS(80).id    = 'NOAA-17 HIRS/3';

KXS(81).value = 228;
KXS(81).id    = 'Japanese Meterological Agency (JMA) profiler U,V';
KXS(82).value = 229;
KXS(82).id    = 'Wind profiler from PILOT (PIBAL) bulletins U,V';
KXS(83).value = 257;
KXS(83).id    = 'NASA/MODIS POES Aqua/Terra IR cloud drift U,V';
KXS(84).value = 258;
KXS(84).id    = 'NASA/MODIS POES Aqua/Terra Water Vapor cloud top U,V';
KXS(85).value = 259;
KXS(85).id    = 'NASA/MODIS POES Aqua/Terra Water Vapor deep layer U,V';

KXS(86).value = 49;
KXS(86).id    = 'Aqua AIRS';
KXS(87).value = 349;
KXS(87).id    = 'Aqua AMSUA';

KXS(88).value = 466;
KXS(88).id    = 'NOAA-16 SBUV/2';


% sensor information:
% ------------------
SENSORS( 1).id  = 'HIRS';
SENSORS( 1).kxs = [14 16 17];
SENSORS( 1).chns = 1:19;
SENSORS( 2).id  = 'MSU';
SENSORS( 2).kxs = [214];
SENSORS( 2).chns = 1:4;
SENSORS( 3).id  = 'SSU';
SENSORS( 3).kxs = [];
SENSORS( 3).chns = [];
SENSORS( 4).id  = 'AMSUA';
SENSORS( 4).kxs = [315 316 349];
SENSORS( 4).chns = 1:15;
SENSORS( 5).id  = 'AMSUB';
SENSORS( 5).kxs = [415 416 417];
SENSORS( 5).chns = 1:5;
SENSORS( 6).id  = 'SSMI';
SENSORS( 6).kxs = [];
SENSORS( 6).chns = [];
SENSORS( 7).id  = 'AIRS';
SENSORS( 7).kxs = [49];
SENSORS( 7).chns = 1:281;
SENSORS( 8).id  = 'METEOSAT';
SENSORS( 8).kxs = [];
SENSORS( 8).chns = [];
SENSORS( 9).id  = 'GOES Sounder';
SENSORS( 9).kxs = [60 62];
SENSORS( 9).chns = 1:18;

% qc history marks:
% ----------------
QCHS(1).value = 0;
QCHS(1).id    = 'none';

% qc exclusion marks:
% ------------------
QCXS(1).value = 0;
QCXS(1).id    = 'none' ;
QCXS(2).value = 1;
QCXS(2).id    = 'passive';
QCXS(3).value = 2;
QCXS(3).id    = 'rejected by GSI';
QCXS(4).value = 3;
QCXS(4).id    = 'very large sigo';
QCXS(5).value = 4;
QCXS(5).id    = 'omf/oma undefined';

% Magnitudes larger than MAXR are replaced by Infs:
% ------------------------------------------------
MAXR = 1e9; 
