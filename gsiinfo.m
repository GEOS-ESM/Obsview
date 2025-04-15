OBSATTRIBUTES.names = {'kx','ks','kt','time','lat','lon','lev',...
    'obs','sigo','omf','oma','qch','qcx','xm'};
OBSATTRIBUTES.descr = {'data source','sounding index','data type',...
    'observation time','latitude','longitude','level/channel',...
    'observed value','observation error',...
    'background departure','analysis departure',...
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
KTS(7).msigo = 10;
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

% Kovach
% Ocean datatype descriptions and units:
% --------------------------------------
KTS(11).value = 101;
KTS(11).id    = 'Sub-surface temperature';
KTS(11).units = 'C';
KTS(11).msigo = 100;
KTS(12).value = 102;
KTS(12).id    = 'Sub-surface salinity';
KTS(12).units = '';
KTS(12).msigo = 100;
KTS(13).value = 104; % 103 in V5
KTS(13).id    = 'Sub-surface zonal velocity';
KTS(13).units = 'm/s';
KTS(13).msigo = 100;
KTS(14).value = 105; % 104 in V5
KTS(14).id    = 'Sub-surface meridional velocity';
KTS(14).units = 'm/s';
KTS(14).msigo = 100;
KTS(15).value = 103; % 105 in V5
KTS(15).id    = 'Sea-surface height anomaly';
KTS(15).units = 'meters';
KTS(15).msigo = 100;
KTS(16).value = 106;
KTS(16).id    = 'Synthetic Salinity';
KTS(16).units = '';
KTS(16).msigo = 100;

KTS(17).value = 17;
KTS(17).id    = 'Rain Rate';
KTS(17).units = 'mm/hr';
KTS(17).msigo = 100;
KTS(18).value = 18;
KTS(18).id    = 'Total Precipitable Water';
KTS(18).units = '';
KTS(18).msigo = 100;
KTS(19).value = 88;
KTS(19).id    = 'Refractivity';
KTS(19).units = 'N';
KTS(19).msigo = 100;
KTS(20).value = 43;
KTS(20).id    = 'Log transformed AOD';
KTS(20).units = 'N';
KTS(20).msigo = 10000;


% pressure-level data ('lev' attribute = pressure):
% -------------------
KTPRS = [4 5 11 21 22 44 88];

% Kovach
% Depth-level data ('lev' attribute = depth);
% ----------------
KTDPH = [101 102 106];

% surface data ('lev' attribute not used):
% ------------
KTSFC = [12 33 39 43 103 104 105];

% radiance data ('lev' attribute = channel number):
% -------------
KTRAD = [40];

% aerosol data ('lev' attribute = channel number):
% -------------
KTAER = [43];

% wind vectors
% ------------
KTWND = {[4 5]};

% data source descriptions  UPDATED 20110503  
% ------------------------------------------
% (from http://www.emc.ncep.noaa.gov/mmb/data_processing/prepbufr.doc/table_2.htm)

KXS(1).value = 301;
KXS(1).id    = 'TERRA MODIS Aerosol (Dark Target Ocean Algorithm)';
KXS(2).value = 302;
KXS(2).id    = 'TERRA MODIS Aerosol (Dark Target Land Algorithm)';
KXS(3).value = 311;
KXS(3).id    = 'AQUA MODIS Aerosol (Dark Target Ocean Algorithm)';
KXS(4).value = 312;
KXS(4).id    = 'AQUA MODIS Aerosol (Dark Target Land Algorithm)';
KXS(5).value = 313;
KXS(5).id    = 'MISR (Multi-angle Imaging SpectroRadiometer)';
KXS(6).value = 314;
KXS(6).id    = 'OMI (Ozone Monitoring Instrument)';
KXS(7).value = 323;
KXS(7).id    = 'AERONET';
KXS(8).value = 310;
KXS(8).id    = 'TERRA MODIS Aerosol (Deep Blue Land Algorithm)';
KXS(9).value = 320;
KXS(9).id    = 'AQUA MODIS Aerosol (Deep Blue Land Algorithm)';


% sensor information:
% ------------------
SENSORS(1).id  = 'HIRS';
SENSORS(1).kxs = [5 6 7 8 9 10 11 12 14 15 16 17 18 19 25];
SENSORS(1).chns = 1:19;


% qc history marks:
% ----------------
QCH_INFLATE = 1.5;
f = 1;
QCHS(1).value = 0;
QCHS(1).id    = 'None';
for k = 2:4
    QCHS(k).value = k-1;
    QCHS(k).id    = ['(' num2str(f,2) ',' num2str(f*QCH_INFLATE,2) '] sigo inflation'];
    f = f*QCH_INFLATE;
end
QCHS(5).value = 4;
QCHS(5).id    = ['(' num2str(f,2) ',' num2str(Inf) '] sigo inflation'];

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
QCXS(5).value = 7;
QCXS(5).id    = 'passive';

% Magnitudes larger than MAXR are replaced by Infs:
% ------------------------------------------------
MAXR = 1e9; 
