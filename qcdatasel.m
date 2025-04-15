function selection = qcdatasel(selid)

% QCDATASEL - Defines data selection

selection = struct('id' ,[],... % identifying string
                   'kt' ,[],... % row vector of data type indices
                   'kx' ,[],... % row vector of data source indices
                   'lev',[],... % row vector of pressure levels
                   'reg',[]);   % region identification string

isel = 0;

isel = isel + 1;
selection(isel).id  = 'All data';

isel = isel + 1;
selection(isel).id  = 'Rawinsonde heights';
selection(isel).kt  = 6;
selection(isel).kx  = 7;

isel = isel + 1;
selection(isel).id  = 'Rawinsonde winds';
selection(isel).kt  = [4 5];
selection(isel).kx  = 7;

isel = isel + 1;
selection(isel).id  = 'Rawinsonde moisture';
selection(isel).kt  = 7;
selection(isel).kx  = 7;

isel = isel + 1;
selection(isel).id  = 'TOVS heights';
selection(isel).kt  = 6;
selection(isel).kx  = [33:56 93:113 125:145 186:227];

isel = isel + 1;
selection(isel).id  = 'TOVS moisture';
selection(isel).kt  = 7;
selection(isel).kx  = [33:56 93:113 125:145 186:227];

isel = isel + 1;
selection(isel).id  = 'SSMI winds';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = [148 149];

isel = isel + 1;
selection(isel).id  = 'QSCAT winds';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = [154];

isel = isel + 1;
selection(isel).id  = 'Surface station sea level pressure';
selection(isel).kt  = [3 6];
selection(isel).kx  = [1 2];

isel = isel + 1;
selection(isel).id  = 'Surface station METAR sea level pressure';
selection(isel).kt  = [3 6];
selection(isel).kx  = [90];

isel = isel + 1;
selection(isel).id  = 'Ship-1 sea level pressure';
selection(isel).kt  = [3 6];
selection(isel).kx  = 3;

isel = isel + 1;
selection(isel).id  = 'Ship-2 sea level pressure';
selection(isel).kt  = [3 6];
selection(isel).kx  = 4;

isel = isel + 1;
selection(isel).id  = 'Environment buoy sea level pressure';
selection(isel).kt  = [3 6];
selection(isel).kx  = 5;

isel = isel + 1;
selection(isel).id  = 'Environment buoy winds';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = 5;

isel = isel + 1;
selection(isel).id  = 'Drifting buoy sea level pressure';
selection(isel).kt  = [3 6];
selection(isel).kx  = 6;

isel = isel + 1;
selection(isel).id  = 'Drifting buoy winds';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = 6;

isel = isel + 1;
selection(isel).id  = 'Ship-1 winds';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = 3;

isel = isel + 1;
selection(isel).id  = 'Ship-2 winds';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = 4;

isel = isel + 1;
selection(isel).id  = 'Dropwinsonde heights';
selection(isel).kt  = 6;
selection(isel).kx  = 11;

isel = isel + 1;
selection(isel).id  = 'Dropwinsonde winds';
selection(isel).kt  = [4 5];
selection(isel).kx  = 11;

isel = isel + 1;
selection(isel).id  = 'Pilot winds';
selection(isel).kt  = [4 5];
selection(isel).kx  = 8;

isel = isel + 1;
selection(isel).id  = 'ACARS';
selection(isel).kt  = [4 5];
selection(isel).kx  = 89;

isel = isel + 1;
selection(isel).id  = 'Aircraft-Air/Sat winds';
selection(isel).kt  = [4 5];
selection(isel).kx  = 14;

isel = isel + 1;
selection(isel).id  = 'Aircraft-Rep winds';
selection(isel).kt  = [4 5];
selection(isel).kx  = 16;

isel = isel + 1;
selection(isel).id  = 'Aircraft-EU winds';
selection(isel).kt  = [4 5];
selection(isel).kx  = 32;

isel = isel + 1;
selection(isel).id  = 'Cloud track winds NESS';
selection(isel).kt  = [4 5];
selection(isel).kx  = 24:25;

isel = isel + 1;
selection(isel).id  = 'Cloud track winds EU';
selection(isel).kt  = [4 5];
selection(isel).kx  = [26 267 268 275:279 286];

isel = isel + 1;
selection(isel).id  = 'Cloud track winds Japanese';
selection(isel).kt  = [4 5];
selection(isel).kx  = [27 269 270];

isel = isel + 1;
selection(isel).id  = 'Cloud track winds US';
selection(isel).kt  = [4 5];
selection(isel).kx  = 119:124;

isel = isel + 1;
selection(isel).id  = 'ERS';
selection(isel).kt  = [1 2 4 5];
selection(isel).kx  = [151 152];

isel = isel + 1;
selection(isel).id  = 'Pseudoheights';
selection(isel).kt  = 6;
selection(isel).kx  = 87;

