function s = obsgroups

% OBSGROUPS Defines observation groups.

% Edit at will. Supported fields are:
%
% s(k).id:   description
% s(k).kx:   data source indices  (enumeration)
% s(k).kt:   data type indices    (enumeration)
% s(k).lev:  levels               (range)
% s(k).lat:  latitudes            (range)
% s(k).lon:  longitudes           (range)
% s(k).qcx:  qc exclusion marks   (enumeration)
% s(k).qch:  qc history marks     (enumeration)
%
% A missing or empty field means no selection on this attribute.

k = 0;

plev = [1000 925 850 700 500 400 300 250 200 150 100 70 50 30 20 10 7 5];

k = k + 1;
s(k).id = 'All data';
s(k).ptstats = false;
s(k).pcoverg = false;

k = k + 1;
s(k).id = 'TERRA MODIS Aerosol (Dark Target Ocean Algorithm)';
s(k).kt = 43;
s(k).kx = 301;
s(k).plev = plev;

k = k + 1;
s(k).id = 'TERRA MODIS Aerosol (Dark Target Land Algorithm)';
s(k).kt = 43;
s(k).kx = 302;
s(k).plev = plev;

k = k + 1;
s(k).id = 'TERRA MODIS Aerosol (Deep Blue Land Algorithm)';
s(k).kt = 43;
s(k).kx = 310;
s(k).plev = plev;

k = k + 1;
s(k).id = 'AQUA MODIS Aerosol (Dark Target Ocean Algorithm)';
s(k).kt = 43;
s(k).kx = 311;
s(k).plev = plev;

k = k + 1;
s(k).id = 'AQUA MODIS Aerosol (Dark Target Land Algorithm)';
s(k).kt = 43;
s(k).kx = 312;
s(k).plev = plev;

k = k + 1;
s(k).id = 'AQUA MODIS Aerosol (Deep Blue Land Algorithm)';
s(k).kt = 43;
s(k).kx = 320;
s(k).plev = plev;

k = k + 1;
s(k).id = 'MISR (Multi-angle Imaging SpectroRadiometer)';
s(k).kt = 43;
s(k).kx = 313;
s(k).plev = plev;

k = k + 1;
s(k).id = 'OMI (Ozone Monitoring Instrument)';
s(k).kt = 43;
s(k).kx = 314;
s(k).plev = plev;

k = k + 1;
s(k).id = 'AERONET';
s(k).kt = 43;
s(k).kx = 323;
s(k).plev = plev;
