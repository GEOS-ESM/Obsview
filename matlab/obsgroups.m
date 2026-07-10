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
s(k).id = 'All temperature data';
s(k).kt = [39 44];
s(k).ptstats = false;
s(k).pcoverg = false;
k = k + 1;
s(k).id = 'Radiosonde temperatures';
s(k).kt = 44;
s(k).kx = 120;
s(k).plev = plev;
k = k + 1;
s(k).id = 'Aircraft temperatures (MDCARS)';
s(k).kt = 44;
s(k).kx = 133;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Aircraft temperatures (SDAR)';
s(k).kt = 44;
s(k).kx = 131;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Aircraft temperatures (AIREP, PIREP)';
s(k).kt = 44;
s(k).kx = 130;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Dropsonde temperatures';
s(k).kt = 44;
s(k).kx = 132;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Surface marine temperatures';
s(k).kt = 44;
s(k).kx = 180;
k = k + 1;
s(k).id = 'Splash level temperatures';
s(k).kt = 44;
s(k).kx = 182;
k = k + 1;
s(k).id = 'Sea-surface temperature data';
s(k).kt = 39;

k = k + 1;
s(k).id = 'All surface pressure data';
s(k).kt = 33;
s(k).ptstats = false;
s(k).pcoverg = false;
k = k + 1;
s(k).id = 'Metar surface pressures';
s(k).kt = 33;
s(k).kx = 187;
k = k + 1;
s(k).id = 'Land station surface pressures';
s(k).kt = 33;
s(k).kx = 181;
k = k + 1;
s(k).id = 'Marine surface pressures';
s(k).kt = 33;
s(k).kx = 180;
k = k + 1;
s(k).id = 'Radiosonde station surface pressures';
s(k).kt = 33;
s(k).kx = 120;
k = k + 1;
s(k).id = 'Splash level surface pressures';
s(k).kt = 33;
s(k).kx = 182;

k = k + 1;
s(k).id = 'All wind data';
s(k).kt = [4 5 12];
s(k).ptstats = false;
s(k).pcoverg = false;
k = k + 1;
s(k).id = 'Radiosonde wind vectors';
s(k).kt = [4 5];
s(k).kx = 220;
s(k).plev = plev;
k = k + 1;
s(k).id = 'Pibal wind vectors';
s(k).kt = [4 5];
s(k).kx = 221;
s(k).plev = plev;
k = k + 1;
s(k).id = 'Wind profiler wind vectors';
s(k).kt = [4 5];
s(k).kx = 223;
s(k).plev = plev;
k = k + 1;
s(k).id = 'NEXRAD wind vectors';
s(k).kt = [4 5];
s(k).kx = 224;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Aircraft wind vectors (MDCARS)';
s(k).kt = [4 5];
s(k).kx = 233;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Aircraft wind vectors (ASDAR)';
s(k).kt = [4 5];
s(k).kx = 231;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Aircraft wind vectors (AIREP, PIREP)';
s(k).kt = [4 5];
s(k).kx = 230;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Dropsonde wind vectors';
s(k).kt = [4 5];
s(k).kx = 232;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Cloud drift wind vectors (METEOSAT-5, METEOSAT-7)';
s(k).kt = [4 5];
s(k).kx = [243 253];
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Cloud drift wind vectors (GMS-5)';
s(k).kt = [4 5];
s(k).kx = [242 252];
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Cloud drift wind vectors (GOES-8 IR, GOES-10 IR)';
s(k).kt = [4 5];
s(k).kx = 245;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Cloud drift wind vectors (GOES-8 WV, GOES-10 WV)';
s(k).kt = [4 5];
s(k).kx = 246;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'PILOT wind vectors';
s(k).kt = [4 5];
s(k).kx = 229;
s(k).plev = plev(plev>=200);
k = k + 1;
s(k).id = 'Surface marine wind vectors';
s(k).kt = [4 5];
s(k).kx = 280;
k = k + 1;
s(k).id = 'QuikSCAT wind vectors';
s(k).kt = [4 5];
s(k).kx = 285;
k = k + 1;
s(k).id = 'ATLAS buoy wind vectors';
s(k).kt = [4 5];
s(k).kx = 282;
k = k + 1;
s(k).id = 'SSM/I wind speeds';
s(k).kt = 12;
s(k).kx = 283;

k = k + 1;
s(k).id = 'All moisture data';
s(k).kt = 11;
s(k).ptstats = false;
s(k).pcoverg = false;
k = k + 1;
s(k).id = 'Radiosonde specific humidities';
s(k).kt = 11;
s(k).kx = 120;
s(k).plev = plev(plev>=300);
k = k + 1;
s(k).id = 'Surface marine specific humidities';
s(k).kt = 11;
s(k).kx = 180;
k = k + 1;
s(k).id = 'Splash level specific humidities';
s(k).kt = 11;
s(k).kx = 182;

k = k + 1;
s(k).id = 'All ozone data';
s(k).kt = [21 22];
s(k).ptstats = false;
s(k).pcoverg = false;
k = k + 1;
s(k).id = 'NOAA-16 SBUV/2 ozone layers';
s(k).kt = 22;
s(k).kx = 466;
s(k).plev = [1013 253 127 63 31 16 7.9 4.0 2.0 1.0 0.50 0.24];
k = k + 1;
s(k).id = 'NOAA-16 SBUV/2 total column ozone';
s(k).kt = 21;
s(k).kx = 466;

k = k + 1;
s(k).id = 'All radiance data';
s(k).kt = 40;
s(k).ptstats = false;
s(k).pcoverg = false;

[SENSORS,KXS] = dconfig('SENSORS','KXS');
for i = 1:length(SENSORS)
   if length(SENSORS(i).chns)<30,
      plev = SENSORS(i).chns;
   else
      [n,p] = hist(SENSORS(i).chns,20);
      plev = round(p);
   end
   for kx = SENSORS(i).kxs,
      k = k + 1;
      s(k).id = [KXS(kx==[KXS.value]).id ' brightness temperatures'];
      s(k).kt = 40;
      s(k).kx = kx;
      s(k).plev = plev;
   end
end

k = k + 1;
s(k).id = 'All synthetic data';
s(k).kx = [111 210];
