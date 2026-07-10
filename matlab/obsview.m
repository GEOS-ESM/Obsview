function obsview(arg,varargin)

% OBSVIEW View observing system.
%
% OBSVIEW without arguments will bring up a file selection dialogue.
% OBSVIEW ODSFILE will open the file ODSFILE.

% 19Oct2001 Dick Dee (dee@dao.gsfc.nasa.gov)
% 23Apr2004 Dick Dee - ECMWF version
% 18Nov2004 Dick Dee - GSI version

warning on

if nargin&&~isodsfile(arg), % must be a callback
    feval(arg,varargin{:})      % simply pass to the callback routine..
    return                      % .. and get out.
end

% must be a user call
% -------------------
if any(findobj('Tag','ObsviewGui')), return; end   % can't have more than one

% create the gui window:
% ---------------------
hGui = figure('Tag','ObsviewGui','Visible','off',...
        'Units','characters','Position',[10 5 80 7.5],...
        'DockControls','off',...
        'IntegerHandle','off','Resize','off','Color',[0.95 0.95 0.8], ...
        'NumberTitle','off','Name','Obsview','MenuBar','none');
    
% set up the menubar:
% ------------------
h = uimenu(hGui,'Label','File');
uimenu(h,'Label','Open ods file',...
    'Callback','obsview openods');
uimenu(h,'Label','Generate html','Separator','on',...
    'Callback','obsview makehtml');
uimenu(h,'Label','Close all figures','Separator','on',...
    'Callback','close(setdiff(findobj(''Type'',''figure''),gcf))');
uimenu(h,'Label','Exit',...
    'Callback','close(gcf)');

h = uimenu(hGui,'Label','Groups');
uimenu(h,'Label','Load new group definitions...',...
    'Callback','obsview opengroups file');
uimenu(h,'Label','Default group definitions',...
    'Callback','obsview opengroups default');
uimenu(h,'Label','Define groups by data source','Separator','on',...
    'Callback','obsview opengroups kx');
uimenu(h,'Label','Define groups by data type',...
    'Callback','obsview opengroups kt');
uimenu(h,'Label','Save group definitions...','Separator','on',...
    'Callback','obsview savegroups');


h = uimenu(hGui,'Label','Coloring');
attr = dconfig('OBSATTRIBUTES'); 
for i = 1:length(attr.names)
    uimenu(h,'Tag',attr.names{i},...
        'Label',['According to ' attr.descr{i}],...
        'Callback',['obsview coloring attr ' attr.names{i}],...
        'Checked','off');
end
uimenu(h,'Separator','on',...
         'Tag','bwr','Label','Blue-white-red',...
    'Callback','obsview coloring cmap bwr')
uimenu(h,'Tag','jet','Label','Blue-cyan-yellow-orange-red',...
    'Callback','obsview coloring cmap jet')
uimenu(h,'Tag','cool','Label','Cyan-magenta',...
    'Callback','obsview coloring cmap cool')
uimenu(h,'Tag','spring','Label','Magenta-yellow',...
    'Callback','obsview coloring cmap spring')
uimenu(h,'Tag','summer','Label','Green-yellow',...
    'Callback','obsview coloring cmap summer')
uimenu(h,'Tag','gray','Label','Grayscale',...
    'Callback','obsview coloring cmap gray')
uimenu(h,'Separator','on',...
         'Tag','keepcolors','Label','Preserve colors when subsetting',...
    'Callback','obsview toggle keepcolors','Checked','on')

h = uimenu(hGui,'Label','Options');
uimenu(h,'Tag','newfig','Label','Always create a new figure when plotting',...
    'Callback','obsview toggle newfig','Checked','on');
uimenu(h,'Tag','qcxclear','Label','Only plot data that passed QC',...
    'Callback','obsview toggle qcxclear','Checked','off');
uimenu(h,'Tag','clegend','Label','Show color legend on plots',...
    'Callback','obsview toggle clegend','Checked','off');
uimenu(h,'Tag','dateline','Label','Center global maps at Int''l Date Line',...
    'Callback','obsview toggle dateline','Checked','off');
s = ver('map'); if ~isempty(s), onoroff = 'on'; else, onoroff = 'off'; end
uimenu(h,'Tag','mapproj','Label','Use map projections for interactive maps',...
    'Callback','obsview toggle mapproj','Enable',onoroff,'Checked',onoroff);

    
% set up the controls:
% -------------------
uicontrol('Tag','Flist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[2 4 76 2], ...
        'Callback','obsview gettimes');
uicontrol('Tag','Tlist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[2 1 22 2]);
uicontrol('Tag','Glist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[26 1 32 2]);
uicontrol('Style','pushbutton','BackgroundColor',[0.95 0.85 0.95], ...
        'Units','characters','Position',[60 1.4 7.8 1.7], ...
        'String','List', ...
        'Callback','obsview showgroup list');
uicontrol('Style','pushbutton','BackgroundColor',[0.8 0.95 0.95], ...
        'Units','characters','Position',[70 1.4 7.8 1.7], ...
        'String','Plot', ...
        'Callback','obsview showgroup plot');

% open an ods file:
% ----------------
if nargin,
    openods(arg)
else,
    openods
end

% open the default group file:
% ---------------------------
opengroups('file',which('obsgroups'))

% set the coloring attr and map
% -------------------------------
coloring('attr','kx')
coloring('cmap','jet')

% make the GUI visible:
% --------------------
set(hGui,'Visible','on','HandleVisibility','callback')

%------------------------------------------------------------------------
function openods(odsfile)

hGui = findobj('Tag','ObsviewGui');
hFlist = findobj(hGui,'Tag','Flist');

if nargin==0,       % get ods file name
    filterspec = '*.ods*';
    path = getappdata(hFlist,'path');
    if any(path),
        filterspec = fullfile(path,filterspec);
    end
    [fname, path] = uigetfile(filterspec, 'Select an ODS file:');
    if fname==0, return, end
    odsfile = [path fname];
end
if ~isodsfile(odsfile), 
    warndlg([odsfile ' is not a readable ods file.'])
    return
end

% store path for next time
% ------------------------
setappdata(hFlist,'path',fileparts(odsfile))

% update gui list of ods files:
% ----------------------------
flist = get(hFlist,'string');
if isempty(deblank(flist{1})), n = 1; else, n = length(flist) + 1; end
flist{n} = odsfile;
set(hFlist,'String',flist,'Value',n)

% update gui list of times on file:
% --------------------------------

gettimes

%------------------------------------------------------------------------
function gettimes

hGui = findobj('Tag','ObsviewGui');
hFlist = findobj(hGui,'Tag','Flist');
flist = get(hFlist,'String'); 
odsfile = flist{get(hFlist,'Value')};  % file name

[fjday,ljday,lhour,nhour,ndata] = getodstimeinfo(odsfile);
dh = 24/nhour;
i = 0;
if sum(ndata>0)>1
    i = i + 1;
    jdays(i) = -1;
    hours(i) = -1;
    tlist{i} = 'All times';
end
j = 1;
for jday = fjday:ljday, 
    for hour = 0:dh:(24-dh),
        if ndata(j)>0,
            i = i + 1;
            jdays(i) = jday; datestr = jdaystr(jday);
            hours(i) = hour; hourstr = [int2str(hour) 'UTC'];
            tlist{i} = [datestr ' ' hourstr];
        end
        j = j + 1;
        if j>length(ndata), break, end
    end
end

hTlist = findobj(hGui,'Tag','Tlist');
set(hTlist,'String',tlist,'Value',min(2,length(tlist)))
setappdata(hTlist,'jdays',jdays)       % store with gui control for convenience
setappdata(hTlist,'hours',hours)

%------------------------------------------------------------------------
function opengroups(varargin)

option = varargin{1};

hGui = findobj('Tag','ObsviewGui');

switch option
    
    case 'default'
        
        groups = obsgroups;

    case 'file'

        if nargin>1,

            grpfile = varargin{2};

        else,       % get group file name from user

            [fname, path] = uigetfile('*.m*', 'Select a data group definition file:');
            if fname==0, return, end    % user pressed cancel
            grpfile = [path fname];

        end

        try,
            [path,fname] = fileparts(grpfile);
            curdir = pwd;
            cd(path)
            groups = feval(fname);
            cd(curdir)
        catch,
            warndlg([grpfile ' is not a proper group definition file.'])
            groups = obsgroups;   % use the default
        end
        
    case {'kx','kt'}

        attr = option;
        
        hFlist = findobj(hGui,'Tag','Flist');
        flist = get(hFlist,'String');
        odsfile = flist{get(hFlist,'Value')};  % file name

        ods = odsload(odsfile,{attr});
        values = unique(ods.(attr));
        for i = 1:length(values), groups(i).(attr) = values(i); end
        
        switch attr
            
            case 'kx'
                
                KXS = dconfig('KXS');
                for i = 1:length(groups) 
                    groups(i).id = KXS(groups(i).kx==[KXS.value]).id; 
                end
                
            case 'kt'
                
                KTS = dconfig('KTS');
                for i = 1:length(groups) 
                    groups(i).id = KTS(groups(i).kt==[KTS.value]).id; 
                end
                
        end

end

% update gui control:
% ------------------
hGlist = findobj(hGui,'Tag','Glist');
for i = 1:length(groups), glist{i} = groups(i).id; end
set(hGlist,'String',glist,'Value',1)
setappdata(hGlist,'groups',groups)

%------------------------------------------------------------------------
function showgroup(action)


hGui   = gcf;
hFlist = findobj(hGui,'Tag','Flist');
hTlist = findobj(hGui,'Tag','Tlist');
hGlist = findobj(hGui,'Tag','Glist');

hwait = waitbar(0,'Please wait...');
set(hGui,'Pointer','watch')

% get the data from file:
% ----------------------
flist = get(hFlist,'String'); 
odsfile = flist{get(hFlist,'Value')};  % file name

if isempty(deblank(odsfile)),
    warndlg('You must first open an ods file.')
    return
end
if ~isodsfile(odsfile), 
    warndlg([odsfile ' is not a readable ods file.'])
    return
end

tlist = get(hTlist,'String'); 
itime = get(hTlist,'Value');           % time index
jdays = getappdata(hTlist,'jdays');
hours = getappdata(hTlist,'hours');

waitbar(0.1,hwait,'Reading observations from file...');

if jdays(itime)<0 
   ods = odsload(odsfile);
else
   ods = odsload(odsfile,jdays(itime),hours(itime));  
end

waitbar(0.3,hwait,'Subsetting observations...');

% subset the data:
% ---------------
igroup = get(hGlist,'Value');
groups = getappdata(hGlist,'groups');

ods = odssubset(ods,groups(igroup));

waitbar(0.5);

if isempty(ods.kt), warndlg('No data in this group.'), end

% do 'quality control' of ods, such as setting NaNs
% -------------------------------------------------
ods = odsclean(ods);

waitbar(0.6,hwait,'Displaying observations...');

% attach subset definitions
% -------------------------
ods.ssid = [tlist{itime} ' ' groups(igroup).id];
if isfield(groups(igroup),'lat'), ods.sslat = groups(igroup).lat; end
if isfield(groups(igroup),'lon'), ods.sslon = groups(igroup).lon; end
if isfield(groups(igroup),'lev'), ods.sslev = groups(igroup).lev; end
if isfield(groups(igroup),'kt' ), ods.sskt  = groups(igroup).kt;  end
if isfield(groups(igroup),'kx' ), ods.sskx  = groups(igroup).kx;  end
if isfield(groups(igroup),'qcx'), ods.ssqcx = groups(igroup).qcx; end
if isfield(groups(igroup),'qch'), ods.ssqch = groups(igroup).qch; end
        
switch action

    case 'plot'

        obsplot(ods)

    case 'list'

        obslist(ods)

end

close(hwait)

set(hGui,'Pointer','arrow')

%------------------------------------------------------------------------
function addgroup

hGui = findobj('Tag','ObsviewGui');    % handle for obsview gui figure
if isempty(hGui), return, end        % obsview gui does not exist

hGlist = findobj(hGui,'Tag','Glist');  % handle for the gui listbox control with group list     
groups = getappdata(hGlist,'groups');  % current list of groups

title = 'Enter a string id for this group:';
prompt = get(findobj(gcbf,'Tag','title2'),'string');   % line 2 of the current plot title
answer = inputdlg(prompt,title);      
if isempty(answer), return, end        % nothing to do

ods = get(gcbf,'Userdata');            % data associated with current plot

k = length(groups)+1;                  % define the new group
groups(k).id = answer{:};          
groups(k).lat = []; if isfield(ods,'sslat'), groups(k).lat = ods.sslat; end
groups(k).lon = []; if isfield(ods,'sslon'), groups(k).lon = ods.sslon; end
groups(k).lev = []; if isfield(ods,'sslev'), groups(k).lev = ods.sslev; end
groups(k).kt  = []; if isfield(ods,'sskt' ), groups(k).kt  = ods.sskt ; end
groups(k).kx  = []; if isfield(ods,'sskx' ), groups(k).kx  = ods.sskx ; end
groups(k).qcx = []; if isfield(ods,'ssqcx'), groups(k).qcx = ods.ssqcx; end
groups(k).qch = []; if isfield(ods,'ssqch'), groups(k).qch = ods.ssqch; end

glist = get(hGlist,'string');          % add new group description to the gui listbox control
glist{k} = groups(k).id;
set(hGlist,'string',glist,'value',k)
setappdata(hGlist,'groups',groups)

%------------------------------------------------------------------------
function savegroups

% get a file name:
% ---------------
[fname,path] = uiputfile('obsgroups.m','Group definition file name');
grpfile = [path fname];

if fname==0, return; end

% attempt to open the file for writing:
% ------------------------------------
fid = fopen(grpfile,'w');
if fid==-1, warnmsg(['Cannot open file ' grpfile ' for writing.']); return, end

% write the preamble:
% ------------------
fprintf(fid,'function s = obsgroups;\n');
fprintf(fid,'\n');
fprintf(fid,'%% This file generated by %s on %s.\n',mfilename,datestr(now));
fprintf(fid,'%%\n');
fprintf(fid,'%% Edit at will. Supported fields are:\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% s(k).id:   description\n');
fprintf(fid,'%% s(k).kx:   data source indices  (enumeration)\n');
fprintf(fid,'%% s(k).kt:   data type indices    (enumeration)\n');   
fprintf(fid,'%% s(k).lev:  levels               (range)\n');       
fprintf(fid,'%% s(k).lat:  latitudes            (range)\n');
fprintf(fid,'%% s(k).lon:  longitudes           (range)\n');
fprintf(fid,'%% s(k).qcx:  qc exclusion marks   (enumeration)\n');
fprintf(fid,'%% s(k).qch:  qc history marks     (enumeration)\n');
fprintf(fid,'%%\n');
fprintf(fid,'%% A missing or empty field means no selection on this attribute.\n');
fprintf(fid,'\n');

% get the group definitions from the gui:
% --------------------------------------
hGui = findobj('Tag','ObsviewGui');    % handle for obsview gui figure
hGlist = findobj(hGui,'Tag','Glist');  % handle for the gui listbox control with group list
groups = getappdata(hGlist,'groups');  % current list of groups

% write the group definitions to file:
% -----------------------------------
fprintf(fid,'k = 0;\n');
for k = 1:length(groups),
    fprintf(fid,'\n');
    fprintf(fid,'k = k + 1;\n');
    fprintf(fid,'s(k).id = ''%s'';\n',groups(k).id);
    for attrs = {'lat','lon','lev','kt','kx','qcx','qch'}, 
        attr = attrs{1};
        if isfield(groups(k),attr), 
            ss = strrep(num2str(groups(k).(attr)),'  ',' ');
            fprintf(fid,['s(k).' attr ' = [' ss '];\n']);
        end
    end
end

% close the file:
% --------------
fclose(fid);

%------------------------------------------------------------------------
function coloring(property,value)

hGui = findobj('Tag','ObsviewGui');
coloring = getappdata(hGui,'Coloring');

if isfield(coloring,property), 
    set(findobj(hGui,'Tag',coloring.(property)),'checked','off'); 
end
coloring.(property) = value;
set(findobj(hGui,'Tag',value),'checked','on')

setappdata(hGui,'Coloring',coloring)

%------------------------------------------------------------------------
function toggle(option)

h = findobj(gcf,'Tag',option);
if strcmp(get(h,'Checked'),'on'), set(h,'Checked','off');
else, set(h,'Checked','on'); end

%------------------------------------------------------------------------
function makehtml

watchon;

hGui = findobj('Tag','ObsviewGui');
hFlist = findobj(hGui,'Tag','Flist');
flist = get(hFlist,'String'); 
odsfile = flist{get(hFlist,'Value')};  % file name

groups = getappdata(findobj(hGui,'Tag','Glist'),'groups');

% close all figures except the gui

close(setdiff(findobj('Type','figure'),hGui))

dirname = mfilename;
if ~exist(dirname,'dir'),
    success = mkdir(dirname);
    if ~success, error(['Cannot create directory ' dirname]); end
else                    % directory already exists
    d = dir(dirname);
    if any(~[d.isdir]), % ... and contains files
        b = questdlg({['Directory ' dirname ' already exists and is not empty. '], ...
                       'Is it OK to delete its contents?'});
        if strcmp(b,'Yes'),
            warning off
            delete([dirname filesep '*'])
            warning on
        else
            warndlg(['Save the contents of directory ' pwd filesep dirname ' and then try again.']) 
            watchoff, return
        end
    end
end

options.groups  = groups;
options.dirname = dirname;
options.browser = true;
options.pcoverg = true;
options.pseries = false;
options.ptstats = true;
options.plotfmt = 'png';
options.plotdpi = '120';

obs2html(mfilename,odsfile,options)

[fjday,ljday,lhour,nhour,nobs] = getodstimeinfo(odsfile);
dh = 24/nhour; hours = 0:dh:(24-dh);
begdate = []; nh = 0;
for jday = fjday:ljday
   for hour = hours
      nh = nh + 1;
      dnum = datenum(jdaystr(jday))+hour/24;
      if nobs(nh)>0
         enddate = [datestr(dnum,'yyyy') datestr(dnum,'mm') datestr(dnum,'dd') sprintf('%02d',hour)];
         if isempty(begdate), begdate = enddate; end
      end
      if jday==ljday && hour==lhour, break; end
   end
end

obs2html(mfilename,begdate,enddate,options)

close(setdiff(findobj('Type','figure'),hGui))

watchoff
