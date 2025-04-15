function obs2html(varargin)

%OBS2HTML Produces html to monitor the observing system.
%
%  Run this program in two steps:
%
%  1) Data processing:
%
%     obs2html(expid,odsfiles[,options])
%
%     - read the data in odsfiles
%     - precompute and store statistics for each analysis
%     - optionally produce data coverage plots
%
%  2) Html creation:
%
%     obs2html(expid,begdate,enddate[,options])
%
%     - compute and store time-averaged statistics
%     - optionally produce time series plots
%     - optionally produce time-averaged statistics plots
%     - create html
%
%  Options:
%
%     options.pcoverg   [logical]  plot coverage map for each group/date
%     options.ptstats   [logical]  plot time-averaged statistics for each group
%     options.pseries   [logical]  plot time series for each group/plev
%
%     options.dirname   [string]   directory for storing results
%     options.grpfile   [string]   .m file containing group definitions
%     options.grpnums   [enumeration]   list of group numbers to be processed
%     options.groups    [structure]   group definition structure 
%                      (NOTE: if options.groups is present then both 
%                       options.grpfile and options.grpnums are ignored)    
%     options.dsynhhs   [enumeration]   list of daily analysis hours to be processed
%     options.tsintvl   [integer]  interval (in days) for time series plots
%
%     options.plotfmt   [string]   graphics format for saved plots
%     options.plotdpi   [string]   plot resolution in dots per inch
%     options.browser   [logical]  open browser to view html
%
%  Last updated 20Jul05 by Dick Dee (dee@gmao.gsfc.nasa.gov)

if nargin==0, help(mfilename); return; end

% optional last argument is options structure
% -------------------------------------------
if isstruct(varargin{nargin}) 
   options = varargin{nargin};
   narg = nargin - 1;
else
   options = [];
   narg = nargin;
end

% default options
% ---------------
if ~isfield(options,'grpfile'), options.grpfile = 'obsgroups'; end
if ~isfield(options,'grpnums'), options.grpnums = []         ; end
if ~isfield(options,'dsynhhs'), options.dsynhhs = [0 3 6 9 12 15 18 21]; end
if ~isfield(options,'tsintvl'), options.tsintvl = 1          ; end
if ~isfield(options,'dirname'), options.dirname = 'obsmon'   ; end

if ~isfield(options,'pcoverg'), options.pcoverg = true       ; end
if ~isfield(options,'ptstats'), options.ptstats = true       ; end
if ~isfield(options,'pseries'), options.pseries = true       ; end

if ~isfield(options,'plotfmt'), options.plotfmt = 'png'      ; end
if ~isfield(options,'plotdpi'), options.plotdpi = '120'      ; end
if ~isfield(options,'browser'), options.browser = true       ; end

switch narg

   case 2 % data processing

      expid    = varargin{1};
      odsfiles = varargin{2};
      
      grpfile = options.grpfile;
      grpnums = options.grpnums;
      dirname = options.dirname;
      
      %warning off MATLAB:MKDIR:DirectoryExists
      
      if isfield(options,'groups') && isstruct(options.groups),
         groups = options.groups;
      else

         disp(['Defining observation groups as in ' which(grpfile)])

         % observation group definitions
         % -----------------------------
         if strcmp(grpfile(end-1:end),'.m'), grpfile(end-1:end) = []; end
         try, groups = feval(grpfile);
	 groups
         catch, error([grpfile ' is not a valid .m file.']), end

      end

      % check against saved definitions, if any
      % ---------------------------------------
      fname = [expid '.' mfilename '.mat'];
      fpath = [dirname filesep fname];
      if exist(fpath,'file')
         disp(['Checking consistency with ' fpath])
         gnew = groups;
         load(fpath,'groups')
         if ~isequal(groups,gnew)
            error('Group definitions differ from previously saved group definitions.')
         end
      else    % save these group definitions
         if ~mkdir(dirname), error(['Cannot create directory ' dirname]); end
         disp(['Saving group definitions in ' fpath])
         save(fpath,'groups')
      end

      % optionally use only a subset of the groups
      % ------------------------------------------
      if any(grpnums)
         if min(grpnums)<1 || max(grpnums)>length(groups), error('Invalid group numbers.'); end
      else
         grpnums = 1:length(groups);
      end

      % process ods files one at a time
      % -------------------------------
      nfs = 0;
      if iscell(odsfiles)       % cell array with file names
         for i = 1:length(odsfiles)
            if isodsfile(odsfiles{i})
               nfs = nfs + 1;
               makedata(odsfiles{i},expid,groups,grpnums,options)
            end
         end
      elseif ischar(odsfiles)   % file name, possibly including '*'
         sdir = dir(odsfiles);
         path = stripfn(odsfiles);
         for i = 1:length(sdir)
            odsfile = [path sdir(i).name];
            if isodsfile(odsfile)
               nfs = nfs + 1;
               makedata(odsfile,expid,groups,grpnums,options)
            end
         end
      else
         help(mfilename)
         error('Check your calling arguments.')
      end

      disp(['Processed ' int2str(nfs) ' ODS files.'])

   case 3 % html creation

      expid   = varargin{1};
      begdate = varargin{2};
      enddate = varargin{3};
      
      grpnums = options.grpnums;
      dsynhhs = options.dsynhhs;
      tsintvl = options.tsintvl;      
      dirname = options.dirname;
      ptstats = options.ptstats;
      pseries = options.pseries;
      plotfmt = options.plotfmt;
      browser = options.browser;
      
      % create time array
      % -----------------
      try, t1 = datenum(begdate,'yyyymmddhh');
      catch, error(['Cannot parse begdate=' begdate]), end
      try, t2 = datenum(enddate,'yyyymmddhh');
      catch, error(['Cannot parse enddate=' enddate]), end
      n = 0; dates = [];
      for dd = floor(t1):tsintvl:ceil(t2)
         for hh = dsynhhs
            t = dd + hh/24;
            if t1<=t && t<=t2
               n = n + 1;
               dates(n).dnum = t;
               dates(n).yyyy = datestr(t,'yyyy');
               dates(n).mm   = datestr(t,'mm'  );
               dates(n).dd   = datestr(t,'dd'  );
               dates(n).hh   = sprintf('%02d',hh);
               dates(n).id   = [dates(n).dd datestr(t,'mmm') dates(n).yyyy ', ' dates(n).hh 'Z'];
            end
         end
      end
      nds = length(dates);
      if nds==0, error(['No data in interval (' begdate ',' enddate ')']), end

      % get group definitions from file
      % -------------------------------
      fname = [expid '.' mfilename '.mat'];
      fpath = [dirname filesep fname];
      disp(['Loading group definitions from ' fpath])
      load(fpath,'groups')

      % optionally use only a subset of the groups
      % ------------------------------------------
      if any(grpnums)
         if min(grpnums)<1 || max(grpnums)>length(groups), error('Invalid group numbers.'); end
      else
         grpnums = 1:length(groups);
      end
      
      % do the work
      % -----------
      if pseries||ptstats
         makestats(dates,expid,groups,grpnums,options)
      end
      
      % create html
      % -----------
      disp('Creating html...')
      disp('    various tables...')
      maketables(dirname)
      disp('    group descriptions...')
      makegrpdoc(groups,grpnums,dirname)
      disp(['    writing ' dirname filesep 'index.html...']);
      makeindex(expid,dates,groups,grpnums,dirname,plotfmt)
      disp('Done.')
      
      if browser,
         stat = web([pwd filesep dirname filesep 'index.html'],'-new');
         if stat~=0, errordlg('Unable to launch browser',mfilename), end
      end

   otherwise

      help(mfilename)
      error('Check your calling arguments.')

end

%-------------------------------------------------------------------

function makedata(odsfile,expid,groups,grpnums,options)

% process an ODS file; produce coverage plots + statistics

dirname = options.dirname;
dsynhhs = options.dsynhhs;
pcoverg = options.pcoverg;
plotfmt = options.plotfmt;
plotdpi = options.plotdpi;

% time info from ods file
% -----------------------
try, [fjday,ljday,lhour,nhour,nobs] = getodstimeinfo(odsfile);
catch, error([odsfile ' is not a valid ODS file.']), end

% store time information for later use
% ------------------------------------
n = 0; nh = 0; dates = [];
for jday = fjday:ljday
   for hour = dsynhhs
      nh = nh + 1;
      dnum = datenum(jdaystr(jday))+hour/24;
      if nobs(nh)>0
         n = n + 1;
         dates(n).nobs = nobs(nh);
         dates(n).jday = jday;
         dates(n).hour = hour;
         dates(n).dnum = dnum;
         dates(n).yyyy = datestr(dnum,'yyyy');
         dates(n).mm   = datestr(dnum,'mm'  );
         dates(n).dd   = datestr(dnum,'dd'  );
         dates(n).hh   = sprintf('%02d',hour);
         dates(n).id   = [dates(n).dd datestr(dnum,'mmm') dates(n).yyyy ', ' dates(n).hh 'Z'];
      end
      if jday==ljday && hour==lhour, break; end
   end
end
nds = length(dates);
if nds==0, return; end % nothing to do

% create directory structure as needed
% ------------------------------------
for k = grpnums
   gdir = [dirname filesep 'G' int2str(k)];
   if ~mkdir(gdir), error(['Cannot create directory ' gdir]); end
   yyyy = 0;
   for n = 1:nds
      if dates(n).yyyy~=yyyy
         yyyy = dates(n).yyyy;
         sdir = [gdir filesep 'Y' yyyy];
         if ~mkdir(sdir), error(['Cannot create directory ' sdir]); end
         mm = 0;
         if dates(n).mm~=mm
            mm = dates(n).mm;
            sdir = [sdir filesep 'M' mm];
            if ~mkdir(sdir), error(['Cannot create directory ' sdir]); end
         end
      end
   end
end

disp(['Processing ' odsfile '...']);

% process the data
% ----------------
for n = 1:nds    % for each synoptic time

   yyyymmddhh = [dates(n).yyyy dates(n).mm dates(n).dd dates(n).hh];

   % get the data from file:
   % ----------------------
   ods = odsload(odsfile,dates(n).jday,dates(n).hour);
   ods.filename = odsfile;

   disp([dates(n).id ': ' int2str(dates(n).nobs) ' observations']);

   % set qcx for large sigo, and other clean-up:
   % ------------------------------------------
   ods = odsclean(ods);

   % subset the data for each group:
   % ------------------------------
   for k = grpnums

      o = odssubset(ods,groups(k));
      o.ssid = [dates(n).id ' ' groups(k).id];
      if isfield(groups(k),'lat'), o.sslat = groups(k).lat; end
      if isfield(groups(k),'lon'), o.sslon = groups(k).lon; end
      if isfield(groups(k),'lev'), o.sslev = groups(k).lev; end
      if isfield(groups(k),'kt' ), o.sskt  = groups(k).kt;  end
      if isfield(groups(k),'kx' ), o.sskx  = groups(k).kx;  end
      if isfield(groups(k),'qcx'), o.ssqcx = groups(k).qcx; end
      if isfield(groups(k),'qch'), o.ssqch = groups(k).qch; end
      ndata = length(o.kt);

      disp([blanks(3) groups(k).id ': ' int2str(ndata) ' observations'])

      if ndata>0, % for this group:

         % directory for storing plots etc.

         sdir = [dirname filesep 'G' int2str(k) filesep 'Y' dates(n).yyyy filesep 'M' dates(n).mm];

         % compute statistics:

         dstats = makets(groups(k));
         dstats = makets(dstats,dates(n).dnum);

         lev = single(dstats.lev);
         if isempty(lev) % single ts for all lev
            dstats = makets(dstats,1,o,true(size(o.lev)));
         else             % one ts for each lev
            % replace o.lev by one of lev
            for i = 1:length(o.lev) % this could be slow
               [y,j] = min(abs(log(lev/o.lev(i))));
               o.lev(i) = lev(j);
            end
            for j = 1:length(lev)
               dstats = makets(dstats,j,o,(o.lev==lev(j)));
            end
         end

         fname = [expid '.' mfilename '.' yyyymmddhh '.mat'];
         save([sdir filesep fname],'dstats')       % save the stats
         
         % optionally create a coverage plot
         
         makeplot = pcoverg;
         if isfield(groups(k),'pcoverg') && ~isempty(groups(k).pcoverg) 
            makeplot = makeplot && groups(k).pcoverg; 
         end
         if makeplot   
            obsplot(o,'obs2html')
            fname = [expid '.obsplot.' yyyymmddhh '.' plotfmt];
            print(gcf,['-d' plotfmt],[sdir filesep fname],['-r' plotdpi])  % save the plot
         end

      end

   end

end

%-------------------------------------------------------------------

function makestats(dates,expid,groups,grpnums,options)

% compute and plot time-averaged statistics; plot time series

dirname = options.dirname;
dsynhhs = options.dsynhhs;
ptstats = options.ptstats;
pseries = options.pseries;
plotfmt = options.plotfmt;
plotdpi = options.plotdpi;

nds = length(dates);
if nds==0, return, end
if isempty(grpnums), return, end

% create time series using precomputed statistics
% -----------------------------------------------
disp('Creating time series...');

for k = grpnums
   ts(k) = makets(groups(k));
end

for n = 1:nds

   yyyymmddhh = [dates(n).yyyy dates(n).mm dates(n).dd dates(n).hh];
   disp(['Reading statistics for ' dates(n).id])

   for k = grpnums

      sdir = [dirname filesep 'G' int2str(k) filesep 'Y' dates(n).yyyy filesep 'M' dates(n).mm];
      fname = [expid '.' mfilename '.' yyyymmddhh '.mat'];
      fpath = [sdir filesep fname];
      if ~exist(fpath,'file'),
         ts(k) = makets(ts(k),dates(n).dnum);        % extend time series with zeros
      else
         load(fpath,'dstats')                    % load precomputed stats
         ts(k) = makets(ts(k),dates(n).dnum,dstats); % extend time series
      end

   end

end

tintstr = [dates(1).yyyy dates(1).mm dates(1).dd dates(1).hh '-' ...
           dates(end).yyyy dates(end).mm dates(end).dd dates(end).hh];
regstr = {'Global','NH','TR','SH'};

for k = grpnums

   disp([groups(k).id '...'])

   disp([blanks(4) '...saving statistics'])

   sdir = [dirname filesep 'G' int2str(k)];
   fname = [expid '.' mfilename '.' tintstr '.mat'];
   tsk = ts(k);
   save([sdir filesep fname],'tsk')
   
   if sum(sum(ts(k).nobs(1,:,:)))>0  % plot only if there are data in this group
      
      makeplot = pseries;
      if isfield(groups(k),'pseries') && ~isempty(groups(k).pseries)
         makeplot = makeplot && groups(k).pseries;
      end
      if makeplot
         disp([blanks(4) '...plotting time series'])
         figure(1); set(1,'Color','w')
         for j = 1:max(1,length(ts(k).lev))
            if sum(ts(k).nobs(1,j,:))>0  % plot only if there are data for this level
               ptsbars(1,ts(k),1,j,expid)  % global statistics only
               fname = [expid '.pseries.L' int2str(j) '.' tintstr '.' plotfmt];
               print(1,['-d' plotfmt],[sdir filesep fname],['-r' plotdpi])
            end
         end

      end

      makeplot = ptstats;
      if isfield(groups(k),'ptstats') && ~isempty(groups(k).ptstats)
         makeplot = makeplot && groups(k).ptstats;
      end
      if makeplot
         disp([blanks(4) '...plotting statistics'])
         figure(1); set(1,'Color','w')
         for ir = 1:4 % for each region
            if sum(sum(ts(k).nobs(ir,:,:)))>0  % plot only if there are data for this region
               ptotals(1,ts(k),ir,expid)
               fname = [expid '.ptstats.' regstr{ir} '.' tintstr '.' plotfmt];
               print(1,['-d' plotfmt],[sdir filesep fname],['-r' plotdpi])
            end
         end

      end

   end

end

%-------------------------------------------------------------------

function ts = makets(varargin)

switch nargin

   case 1  % initialize time series structure

      group = varargin{1};

      ts.id = group.id;
      ts.reg = {'Global','N. Hemisphere','Tropics','S. Hemisphere'};
      if isfield(group,'kt') && ~isempty(group.kt),
         [KTPRS,KTRAD] = dconfig('KTPRS','KTRAD');
         ts.isprs = all(ismember(group.kt,KTPRS));
         ts.israd = all(ismember(group.kt,KTRAD));
      else
         ts.isprs = false;
         ts.israd = false;
      end
      if isfield(group,'plev') && ~isempty(group.plev), ts.lev = group.plev;
      else, ts.lev = []; end
      ts.time  = [];
      ts.nobs  = [];
      ts.nana  = [];
      ts.npsv  = [];
      ts.somf  = [];
      ts.somf2 = [];
      ts.soma  = [];
      ts.soma2 = [];
      ts.somfa = [];
      ts.sigo  = [];
      ts.sigo2 = [];
      ts.Jobkg = [];
      ts.Joana = [];

   case 2  % extend time series with zeros

      ts   = varargin{1};
      time = varargin{2};

      nr = length(ts.reg);
      nl = max(1,length(ts.lev));
      nt = length(ts.time);
      ts.time (          nt+1) = time;
      ts.nobs (1:nr,1:nl,nt+1) = 0;
      ts.nana (1:nr,1:nl,nt+1) = 0;
      ts.npsv (1:nr,1:nl,nt+1) = 0;
      ts.somf (1:nr,1:nl,nt+1) = 0;
      ts.somf2(1:nr,1:nl,nt+1) = 0;
      ts.soma (1:nr,1:nl,nt+1) = 0;
      ts.soma2(1:nr,1:nl,nt+1) = 0;
      ts.somfa(1:nr,1:nl,nt+1) = 0;
      ts.sigo (1:nr,1:nl,nt+1) = 0;
      ts.sigo2(1:nr,1:nl,nt+1) = 0;
      ts.Jobkg(1:nr,1:nl,nt+1) = 0;
      ts.Joana(1:nr,1:nl,nt+1) = 0;

   case 3  % extend time series with precomputed stats

      ts    = varargin{1};
      time  = varargin{2};
      s     = varargin{3};

      nr = length(ts.reg);
      nl = max(1,length(ts.lev));
      nt = length(ts.time);
      ts.time (          nt+1) = time;
      ts.nobs (1:nr,1:nl,nt+1) = s.nobs (1:nr,1:nl,1);
      ts.nana (1:nr,1:nl,nt+1) = s.nana (1:nr,1:nl,1);
      ts.npsv (1:nr,1:nl,nt+1) = s.npsv (1:nr,1:nl,1);
      ts.somf (1:nr,1:nl,nt+1) = s.somf (1:nr,1:nl,1);
      ts.somf2(1:nr,1:nl,nt+1) = s.somf2(1:nr,1:nl,1);
      ts.soma (1:nr,1:nl,nt+1) = s.soma (1:nr,1:nl,1);
      ts.soma2(1:nr,1:nl,nt+1) = s.soma2(1:nr,1:nl,1);
      ts.somfa(1:nr,1:nl,nt+1) = s.somfa(1:nr,1:nl,1);
      ts.sigo (1:nr,1:nl,nt+1) = s.sigo (1:nr,1:nl,1);
      ts.sigo2(1:nr,1:nl,nt+1) = s.sigo2(1:nr,1:nl,1);
      ts.Jobkg(1:nr,1:nl,nt+1) = s.Jobkg(1:nr,1:nl,1);
      ts.Joana(1:nr,1:nl,nt+1) = s.Joana(1:nr,1:nl,1);

   case 4  % compute j^th time series entries from ods

      ts  = varargin{1};
      j   = varargin{2};  % level index
      o   = varargin{3};  % ods
      msk = varargin{4};  % mask on ods for this level
      
      nr = length(ts.reg);
      nt = length(ts.time);
      for ir = 1:length(ts.reg)     % stats by region
         mskr = msk;
         if     ir==2, mskr = mskr & (o.lat>20);
         elseif ir==3, mskr = mskr & (-20<=o.lat & o.lat<20);
         elseif ir==4, mskr = mskr & (o.lat<-20); 
         end
         ts.nobs(ir,j,nt) = sum(mskr); 
         ana = (o.qcx==0 & mskr); 
         psv = ((o.qcx==1 | o.qcx==7) & mskr); 
         ts.nana(ir,j,nt) = sum(ana);
         ts.npsv(ir,j,nt) = sum(psv);
         i = ana|psv;
         if any(i)
            ts.somf (ir,j,nt) = sum( o.omf(i));
            ts.somf2(ir,j,nt) = sum((o.omf(i)).^2);
            ts.soma (ir,j,nt) = sum( o.oma(i));
            ts.soma2(ir,j,nt) = sum((o.oma(i)).^2);
            ts.somfa(ir,j,nt) = sum( o.omf(i).* o.oma(i) );
         end
         i = ana;
         if any(i)
            ts.sigo (ir,j,nt) = sum( o.sigo(i));
            ts.sigo2(ir,j,nt) = sum((o.sigo(i)).^2);
            ts.Jobkg(ir,j,nt) = sum((o.omf(i)./o.sigo(i)).^2);
            ts.Joana(ir,j,nt) = sum((o.oma(i)./o.sigo(i)).^2);
         end
      end

end

%-------------------------------------------------------------------

function makeindex(expid,dates,groups,grpnums,hdir,plotfmt)

% MAKEINDEX - Produces index.html for obs2html

nds = length(dates);
if nds==0, return; end
tintstr = [dates(1).yyyy dates(1).mm dates(1).dd dates(1).hh '-' ...
           dates(end).yyyy dates(end).mm dates(end).dd dates(end).hh];
regstr = {'Global','NH','TR','SH'};

fid = fopen([hdir filesep 'index.html'],'w');

fprintf(fid,'%s\n','<html>');
fprintf(fid,'%s\n',['<head><title>' expid ' Observing System Summary</title></head>']);
fprintf(fid,'%s\n','<body text="#000000" bgcolor="#ffffe3" link="#0000ff" vlink="#800080" alink="#ff00ff">');
fprintf(fid,'%s\n','<script language="javascript">');
fprintf(fid,'%s\n','<!--');
fprintf(fid,'%s\n','function leapto(form)');
fprintf(fid,'%s\n','{var pname = form.pname.selectedIndex;');
fprintf(fid,'%s\n','window.location = (form.pname.options[pname].value);}');
fprintf(fid,'%s\n','//-->');
fprintf(fid,'%s\n','</script>');

fprintf(fid,'%s\n',['<center><h2>Observing System Summary - ' expid '</h2></center>']);
if length(dates)==1,
   str = dates(1).id;
else
   str = [dates(1).id ' - ' dates(end).id];
end
fprintf(fid,'%s\n',['<center><b><font color=blue>' str '</font></b></center>']);

fprintf(fid,'%s\n','<p><table cellspacing=10>');

fprintf(fid,'%s\n','<tr>');
fprintf(fid,'%s\n','<th align=right></th>');
fprintf(fid,'%s\n','<th align=left></th>');
fprintf(fid,'%s\n','<th colspan=4><font color=red>Statistics</font></th>');
fprintf(fid,'%s\n','<th><font color=red>Time series</font></th>');
fprintf(fid,'%s\n','<th><font color=red>Spatial coverage</font></th>');
fprintf(fid,'%s\n','</tr>');

[KTPRS,KTRAD] = dconfig('KTPRS','KTRAD');

for k=grpnums  % for each group

   gdir = ['G' int2str(k)];

   fprintf(fid,'%s\n','<tr align=right>'); % start a new row
   
   % group description
   
   fprintf(fid,'%s\n',['<td align=right><a href="index.html" ONCLICK=window.open("G' int2str(k) ...
      '.html","","width=500,height=400,resizable,scrollbars")><b>G' int2str(k) ...
      '</b></a></td>']);
   
   if strncmp(groups(k).id,'All ',4)
      fprintf(fid,'%s\n',['<td align=left><b><font color=red>' groups(k).id '</font></b></td>']);
   else
      fprintf(fid,'%s\n',['<td align=left><font color=blue>' groups(k).id '</font></td>']);
   end

   % time-averaged statistics

   for ir = 1:4
      fpath = [gdir filesep expid '.ptstats.' regstr{ir} '.' tintstr '.' plotfmt];
      if exist([hdir filesep fpath],'file')
         fprintf(fid,'%s\n',['<td align=center><a href="' fpath '"><b>' regstr{ir} '</b></a></td>']);
      else
         fprintf(fid,'%s\n','<td><hr></td>');
      end
   end

   % time series by level

   np = 0;
   if ~isfield(groups(k),'plev') || isempty(groups(k).plev)
      j = 1;
      fpath = [gdir filesep expid '.pseries.L' int2str(j) '.' tintstr '.' plotfmt];
      if exist([hdir filesep fpath],'file'),
         fprintf(fid,'%s\n',['<td align=center><a href="' fpath '"><b>Global</b></a></td>']);
         np = np + 1;
      end
   else
      kt = [];
      if isfield(groups(k),'kt'), kt = groups(k).kt; end
      if any(kt) && all(ismember(kt,KTPRS)),     selstr = 'Binned by level:  ';
      elseif any(kt) && all(ismember(kt,KTRAD)), selstr = 'Binned by channel:';
      else selstr = 'Binned by ''lev'':'; end
      nl = length(groups(k).plev);
      for j=nl:-1:1
         fpath = [gdir filesep expid '.pseries.L' int2str(j) '.' tintstr '.' plotfmt];
         if exist([hdir filesep fpath],'file'),
            if np==0
               fprintf(fid,'%s\n','<td><form><select name="pname" onChange="leapto(this.form)">');
               fprintf(fid,'%s\n',['<option selected value="index.html">' selstr '</option>']);
            end
            fprintf(fid,'%s\n',['<option value="' fpath '">' num2str(groups(k).plev(j)) '</option>']);
            np = np + 1;
         end
      end
      if np>0,
         fprintf(fid,'%s\n','</select></form></td>');
      end
   end
   if np==0
      fprintf(fid,'%s\n','<td><hr></td>');
   end

   % coverage by date

   np = 0;
   for n=1:length(dates)
      sdir = [gdir filesep 'Y' dates(n).yyyy filesep 'M' dates(n).mm];
      yyyymmddhh = [dates(n).yyyy dates(n).mm dates(n).dd dates(n).hh];
      fpath = [sdir filesep expid '.obsplot.' yyyymmddhh '.' plotfmt];
      if exist([hdir filesep fpath],'file'),
         if np==0
            fprintf(fid,'%s\n','<td><form><select name="pname" onChange="leapto(this.form)">');
            fprintf(fid,'%s\n','<option selected value="index.html">Analysis time:</option>');
         end
         fprintf(fid,'%s\n',['<option value="' fpath '">' dates(n).id '</option>']);
         np = np + 1;
      end
   end
   if np>0
      fprintf(fid,'%s\n','</select></form></td>');
   end
   if np==0
      fprintf(fid,'%s\n','<td><hr></td>');
   end

   fprintf(fid,'%s\n','</tr>'); % end of this row

end

fprintf(fid,'%s\n','</td></table>');

fprintf(fid,'%s\n','<p><hr>');

fprintf(fid,'%s\n','<ul>');
fprintf(fid,'%s\n','<li><a href="index.html" ONCLICK=window.open("kxtable.html","","width=800,height=500,resizable,scrollbars")>Data Sources (kx)</a></li>');
fprintf(fid,'%s\n','<li><a href="index.html" ONCLICK=window.open("kttable.html","","width=400,height=500,resizable,scrollbars")>Data Types (kt)</a></li>');
fprintf(fid,'%s\n','<li><a href="index.html" ONCLICK=window.open("qcxtable.html","","width=400,height=500,resizable,scrollbars")>QC Exclusion Marks (qcx)</a></li>');
fprintf(fid,'%s\n','<li><a href="index.html" ONCLICK=window.open("qchtable.html","","width=400,height=500,resizable,scrollbars")>QC History Marks (qch)</a></li>');
fprintf(fid,'%s\n','</ul>');

fprintf(fid,'%s\n','<p><hr>');

fprintf(fid,'%s\n',['<h4><i>Generated by OBS2HTML on ' datestr(now,0) '.</h4>']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

%-------------------------------------------------------------------

function maketables(dirname)

% This program creates:
%
%      kxtable.html
%      kttable.html
%      qchtable.html
%      qcxtable.html

[KXS,KTS,QCHS,QCXS] = dconfig('KXS','KTS','QCHS','QCXS');

% create kxtable.html
% -------------------

fid = fopen([dirname filesep 'kxtable.html'],'w');

fprintf(fid,'%s','<html><head><title>Data Source Index Table</title></head>');
fprintf(fid,'%s','<body text="#000000" bgcolor="#ffffcc" link="#0000ff" vlink="#006030" alink="#ff0000">');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><font color="0000ff"><b>Data Source Index Table</b></font></caption>');
[w,i] = sort([KXS.value]);
for k = i,
   fprintf(fid,['<tr><td>%d</td><td>' KXS(k).id '</td>\n'],KXS(k).value);
end
fprintf(fid,'%s','</center></table>');


fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

% create kttable.html
% -------------------

fid = fopen([dirname filesep 'kttable.html'],'w');

fprintf(fid,'%s','<html><head><title>Data Type Index Table</title></head>');
fprintf(fid,'%s','<body text="#000000" bgcolor="#ffffcc" link="#0000ff" vlink="#006030" alink="#ff0000">');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><font color="0000ff"><b>Data Type Index Table</b></font></caption>');
[w,i] = sort([KTS.value]);
for k = i,
   units = KTS(k).units;
   if strcmp(units,'%'), units = '%%'; end
   fprintf(fid,['<tr><td>%d</td><td>' KTS(k).id ' [' units ']</td>\n'],KTS(k).value);
end
fprintf(fid,'%s','</center></table>');

fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

% create qchtable.html
% -------------------

fid = fopen([dirname filesep 'qchtable.html'],'w');

fprintf(fid,'%s','<html><head><title>Quality Control History Marks</title></head>');
fprintf(fid,'%s','<body text="#000000" bgcolor="#ffffcc" link="#0000ff" vlink="#006030" alink="#ff0000">');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><font color="0000ff"><b>Quality Control History Marks</b></font></caption>');

[w,i] = sort([QCHS.value]);
for k = i,
   fprintf(fid,['<tr><td>%d</td><td>' QCHS(k).id '</td>\n'],QCHS(k).value);
end
fprintf(fid,'%s','</center></table>');

fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

% create qcxtable.html
% -------------------

fid = fopen([dirname filesep 'qcxtable.html'],'w');

fprintf(fid,'%s','<html><head><title>Quality Control Exclusion Marks</title></head>');
fprintf(fid,'%s','<body text="#000000" bgcolor="#ffffcc" link="#0000ff" vlink="#006030" alink="#ff0000">');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><font color="0000ff"><b>Quality Control Exclusion Marks</b></font></caption>');

[x,i] = sort([QCXS.value]);
for k = i,
   fprintf(fid,['<tr><td>%d</td><td>' QCXS(k).id '</td>\n'],QCXS(k).value);
end
fprintf(fid,'%s','</center></table>');

fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

%-------------------------------------------------------------------

function makegrpdoc(groups,grpnums,dirname)

% create an html file for each group, describing the selection
% criteria for that group

[KTS,KXS,QCHS,QCXS] = dconfig('KTS','KXS','QCHS','QCXS');

for k = grpnums,

   group = groups(k);

   fid = fopen([dirname filesep 'G' int2str(k) '.html'],'w'); % create an html file for this group

   fprintf(fid,'%s\n',['<html><head><title>' group.id '</title></head>']);
   fprintf(fid,'%s\n','<body text="#000000" bgcolor="#ffffcc" link="#0000ff" vlink="#006030" alink="#ff0000">');

   fprintf(fid,'%s\n',['<h2><font color=blue>' group.id '</font> are all observations for which</h2>']);
   AND = '';

   if isfield(group,'kx')&&~isempty(group.kx),
      fprintf(fid,'%s\n',['<h4>' AND 'kx <font color=blue>(data source)</font> is one of:</h4>']);
      fprintf(fid,'%s\n','<table>');
      for kx = group.kx,
         fprintf(fid,'<tr><td align=right>%i</td><td>%s</td></tr>\n',kx,KXS(kx==[KXS.value]).id);
      end
      fprintf(fid,'%s\n','</table>');
      AND = '<font color=red>AND </font>';
   end

   if isfield(group,'kt')&&~isempty(group.kt),
      fprintf(fid,'%s\n',['<h4>' AND 'kt <font color=blue>(data type)</font> is one of:</h4>']);
      fprintf(fid,'%s\n','<table>');
      for kt = group.kt,
         fprintf(fid,'<tr><td align=right>%i</td><td>%s</td></tr>\n',kt,KTS(kt==[KTS.value]).id);
      end
      fprintf(fid,'%s\n','</table>');
      AND = '<font color=red>AND </font>';
   end

   if isfield(group,'lev')&&~isempty(group.lev),
      lev = group.lev;
      if diff(lev),
         str = ['in the range [' num2str(lev(1)) 'hPa, ' num2str(lev(2)) 'hPa]'];
      else
         str = [num2str(lev(1)) 'hPa'];
      end
      fprintf(fid,'%s\n',['<h4>' AND 'lev <font color=blue>(data level)</font> is ' str '</h4>']);
      AND = '<font color=red>AND </font>';
   end

   if isfield(group,'lat')&&~isempty(group.lat),
      lat = group.lat;
      if diff(lat),
         str = ['in the range [' num2str(lat(1)) 'deg., ' num2str(lat(2)) 'deg.]'];
      else
         str = [num2str(lat(1)) 'deg.'];
      end
      fprintf(fid,'%s\n',['<h4>' AND 'lat <font color=blue>(latitude)</font> is ' str ' </h4>']);
      AND = '<font color=red>AND </font>';
   end

   if isfield(group,'lon')&&~isempty(group.lon),
      lon = group.lon;
      lon = -180 + mod(lon + 180,360);     % lon in [-180,180)
      if diff(lon),
         if lon(1)>=lon(2), lon(2) = lon(2) + 360; end
         str = ['in the range [' num2str(lon(1)) 'deg., ' num2str(lon(2)) 'deg.]'];
      else
         str = [num2str(lon(1)) 'deg.'];
      end
      fprintf(fid,'%s\n',['<h4>' AND 'lon <font color=blue>(longitude)</font> is ' str ' </h4>']);
      AND = '<font color=red>AND </font>';
   end

   if isfield(group,'qch')&&~isempty(group.qch),
      fprintf(fid,'%s\n',['<h4>' AND 'qchist <font color=blue>(quality control history mark)</font> is one of:</h4>']);
      fprintf(fid,'%s\n','<table>');
      for qch = group.qch,
         fprintf(fid,'<tr><td align=right>%i</td><td>%s</td></tr>\n',qch,QCHS(qch==[QCHS.value]).id);
      end
      fprintf(fid,'%s\n','</table>');
      AND = '<font color=red>AND </font>';
   end

   if isfield(group,'qcx')&&~isempty(group.qcx),
      fprintf(fid,'%s\n',['<h4>' AND 'qcexcl <font color=blue>(quality control exclusion mark)</font> is one of:</h4>']);
      fprintf(fid,'%s\n','<table>');
      for qcx = group.qcx,
         fprintf(fid,'<tr><td align=right>%i</td><td>%s</td></tr>\n',qcx,QCXS(qcx==[QCXS.value]).id);
      end
      fprintf(fid,'%s\n','</table>');
      AND = '<font color=red>AND </font>';
   end

   fprintf(fid,'%s\n','<p><hr>');
   fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
   fprintf(fid,'%s\n','</body></html>');

   fclose(fid);

end

