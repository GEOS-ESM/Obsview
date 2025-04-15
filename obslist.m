function obslist(arg,varargin)

% OBSLIST List observations in a scrollable window.
%
% OBSLIST(ODS) produces a scrollable window with a list of data attributes
%     for all observations contained in the ods structure ODS.
%
% OBSLIST(ODSFILE) loads the data from the ODS file ODSFILE.
%
% OBSLIST without arguments brings up a file selection gui.
%
%     Click on any item in the list to bring up another window with more 
%     descriptive information.  Right-click on the list for various options
%     and actions that can be performed.
%
%     Items can be selected by clicking on the buttons at the top of the
%     list, or by clicking on the sequence numbers in the first column.

% 26Aug2002 Dick Dee (dee@dao.gsfc.nasa.gov)
% 23Apr2004 Dick Dee - ECMWF version
% 18Nov2004 Dick Dee - GSI version


if nargin==0,
    [fname, path] = uigetfile('*.ods*', 'Select an ODS file:');
    if fname==0, return, end
    arg = [path fname];
end

if ischar(arg), % either an ODS file or a callback
    if ~isodsfile(arg),
        feval(arg,varargin{:})      % simply pass to the callback routine..
        return                      % .. and get out.
    else
        ods = getods(arg);
        ods = convertods(ods);
    end
else            % otherwise this will set up a new list window
    ods = arg;  % first arg must be an ods structure
end

% check if we have data to plot

nobs = length(ods.kt);
if nobs==0, return, end

attr = dconfig('OBSATTRIBUTES'); 

% replace empty attributes by NaNs:

nans = NaN*ones(nobs,1);
for i = 1:length(attr.names),
    a = attr.names{i};
    if ~isfield(ods,a), ods.(a) = nans; 
    elseif isempty(ods.(a)), ods.(a) = nans; end
end

% attach selection flag to the observations

ods.select = false(nobs,1); 

% create the figure window:

units = get(0,'Units'); 
set(0,'Units','character'); 
scr = get(0,'ScreenSize'); scrh = scr(4); % screen height in character units
set(0,'Units',units)

nx = length(int2str(nobs));               % width needed for obs numbers
ns = 3;                                   % width needed for slider
dx = 4 + [nx 3 7 3 7 7 7 7 7 7 7 7 3 3 7 ns]; % column widths 
xtabs = cumsum(dx);                       % right edges of columns
nv = sum(dx(2:end-1));                    % width of columns with attribute values
nl = min([40 scrh-24]);                   % nominal number of lines per page
nl = max([2 min([nl nobs])]);             % actual number of lines per page

hFig = figure('Units','character','Position',[20 10 sum(dx) nl+6], ...
    'DockControls','off',...
    'PaperOrientation','landscape','PaperPosition',[0 0.5 11 8],...
    'MenuBar','none','Resize','off','Userdata',ods, ...
    'Color',[0 0 0],'Name',[' ' int2str(nobs) ' observations'], ...
    'DefaultAxesUnits','character',...
    'DefaultUIcontrolUnits','character',...
    'CloseRequestFcn','delete([gcf; findobj(''Tag'',''obsinfo'')])');

% create a button for each attribute:
 
for i = 1:length(attr.names),
    uicontrol('Parent',hFig,'Style','pushbutton',...
        'Position',[xtabs(i) nl+4 dx(i+1) 1.2],...
        'String',attr.names{i},'TooltipString',attr.descr{i},...
        'Callback',['obslist select attr ' attr.names{i}]);
end

% create the axes for the line numbers:

axes('Parent',hFig,'Units','character','Position',[1 2 nx+2 nl+1], ...
    'DefaultTextFontName','default', ...
    'DefaultTextHorizontalAlignment','right',...
    'XLim',[0 1],'YLim',[1 nl],'YDir','reverse','Visible','off', ...
    'Tag','numaxes');

% create the axes for the attribute values:

hListAxes = axes('Parent',hFig,'Position',[nx+2 2 nv nl+1], ...
    'DefaultTextFontName','default', ...
    'DefaultTextHorizontalAlignment','right',...
    'XLim',[1 nv],'YLim',[1 nl],'YDir','reverse','Visible','off', ...
    'ButtonDownFcn','obslist clickonlist','Tag','listaxes');

% create a slider control

if nobs>nl,      
      uicontrol('Parent',hFig,'Style','slider',...
            'Units','characters','Position',[xtabs(end-1)+2 1.75 3.5 nl+2],...
            'ForegroundColor',0.4*[1 1 1],...
            'Value',nobs-39,'Min',1,'Max',nobs-39,'SliderStep',[1 30]/(nobs-1),...
            'Callback','obslist dsplay','Tag','slider');
end

% create menubar

mkmenubar(hFig)

% write the first set of text lines:

figure(hFig)
dsplay(1)
 
%-------------------------------------------------------------------------
function dsplay(n)

% n indicates index of first observation on the page

if nargin==0,
    [h,hFig] = gcbo;      % handles to slider control (h) and figure (hFig)
    n = round(get(h,'Max')-get(h,'Value')+1);  
else,
    hFig = gcf;
    h = findobj(hFig,'Tag','slider');
    set(h,'Value',get(h,'Max')-n+1);
end
hListAxes = findobj(hFig,'Tag','listaxes');  % handle to axes containing obs list
hNumAxes  = findobj(hFig,'Tag','numaxes');   % handle to axes containing obs numbers

ods = get(hFig,'Userdata');
nobs = length(ods.kt);
iobs = n:min([nobs n+39]);

delete(findobj(hListAxes,'Tag','obstext'))  % rewrite text lines
delete(findobj(hNumAxes, 'Tag','obsnum'))   
delete(findobj(hListAxes,'Tag','obsselect')) 

% create text objects for all observations:

dx = 4 + [3 7 3 7 7 7 7 7 7 7 7 3 3 7]; % column widths 
x  = 2 + cumsum(dx);                    % right edges of columns

for is = 1:length(iobs),
    i = iobs(is);
    if ods.select(i),  % highlight selected observations
        patch(0.5+[2 x(end) x(end) 2],is-0.25+0.5*[1 1 -1 -1],0.3*[1 1 1],...
            'Parent',hListAxes,'Clipping','off','Tag','obsselect')
    end
    text(1,is,int2str(i),...
        'Parent',hNumAxes,'Color','w',...
        'ButtonDownFcn','obslist clickonnum','Tag','obsnum');
    if isfield(ods,'cidx'),
        colr = ods.cinfo(ods.cidx(i)).rgb;
    else
        colr = [0.8125 1.0000 0.1875];
    end
    text(x( 1),is,sprintf(  '%3d',ods.kx(i)   ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 2),is,sprintf(  '%7d',ods.ks(i)   ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 3),is,sprintf(  '%3d',ods.kt(i)   ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 4),is,sprintf(  '%6d',ods.time(i) ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 5),is,sprintf('%7.2f',ods.lat(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 6),is,sprintf('%7.2f',ods.lon(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 7),is,sprintf('%7.2f',ods.lev(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 8),is,sprintf('%7.2f',ods.obs(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x( 9),is,sprintf('%6.2f',ods.sigo(i) ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x(10),is,sprintf('%6.2f',ods.omf(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x(11),is,sprintf('%6.2f',ods.oma(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x(12),is,sprintf(  '%3d',ods.qch(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x(13),is,sprintf(  '%3d',ods.qcx(i)  ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
    text(x(14),is,sprintf('%7.2f',ods.xm(i)   ),'Parent',hListAxes,'Color',colr,'Tag','obstext');
end

% store indices of displayed data with the obslist axes

set(hListAxes,'Userdata',iobs)

% transfer mouse action

set(allchild(hListAxes),'UIContextMenu',get(hListAxes,'UIContextMenu'));
set(allchild(hListAxes),'ButtonDownFcn',get(hListAxes,'ButtonDownFcn'));

%-------------------------------------------------------------------------
function select(action,varargin)

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)
ods = get(hFig,'Userdata');   % get ods stored with the list

switch action

   case 'all'

      ods.select = true(size(ods.kt));

   case 'none'

      ods.select = false(size(ods.kt));

   case 'attr'

      attr = varargin{1};

      ss = ods.select;            % sub-select..
      if ~any(ss), ss = ~ss; end  % unless nothing has been selected yet

      ods.select = false(size(ss));

      title = 'Select data attributes';

      [OBSATTRIBUTES,KTS,KXS,QCHS,QCXS] = ...
         dconfig('OBSATTRIBUTES','KTS','KXS','QCHS','QCXS');
      attrs = OBSATTRIBUTES.names;
      descr = OBSATTRIBUTES.descr;
      discr = OBSATTRIBUTES.discr;

      ia = find(strcmp(attr,attrs));
      if any(ia)

         if discr(ia) % discrete attribute values

            values = unique(ods.(attr)(ss)); n = length(values);
            prompt = ['Select one or more ' descr{ia} ':'];
            switch attr
               case 'kt' , for i=1:n, listr{i} = KTS([KTS.value]==values(i)).id; end
               case 'kx' , for i=1:n, listr{i} = KXS([KXS.value]==values(i)).id; end
               case 'qch', for i=1:n, listr{i} = QCHS([QCHS.value]==values(i)).id; end
               case 'qcx', for i=1:n, listr{i} = QCXS([QCXS.value]==values(i)).id; end
               otherwise, listr = int2str(values);
            end
            i = listdlg('Name',title,'PromptString',prompt,...
               'ListSize',[240 120],'ListString',listr);
            if isempty(i), return; end

            for v = values(i)', ods.select = ods.select | (ods.(attr)==v); end

         else        % continuous attribute values

            minx = min(ods.(attr)(ss));
            maxx = max(ods.(attr)(ss));
            prompt  = {['Minimum ' descr{ia} ':'],['Maximum ' descr{ia} ':']};
            def     = {num2str(minx),num2str(maxx)};
            answer  = inputdlg(prompt,title,1,def);
            if isempty(answer), return; end
            minx    = eval(answer{1});
            maxx    = eval(answer{2});

            ods.select = (minx<=ods.(attr) & ods.(attr)<=maxx);

         end

      end

      ods.select = ods.select & ss;

   case 'list'

      switch varargin{1}

         case 'selected'

            obslist(odssubset(ods,ods.select))
            return

         case 'unselected'

            obslist(odssubset(ods,~ods.select))
            return

      end
end

nobs = length(ods.select);
nsel = sum(ods.select);
set(hFig,'Userdata',ods,...
   'Name',[' ' int2str(nobs) ' observations (' int2str(nsel) ' selected)'])

is = find(ods.select);
if nsel, n = is(1); else, n = 1; end
dsplay(n)

%-------------------------------------------------------------------------
function clickonlist

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)

% in case of right-click, do nothing:

if strcmp(get(hFig,'SelectionType'),'alt'), return, end

% get ods stored with the list:

ods = get(hFig,'Userdata');

% get mouse click position and index of the corresponding obs:

hAxes = findobj(hFig,'Tag','listaxes');  % handle to axes containing obs list
v = get(hAxes,'CurrentPoint');
iy = round(v(1,2));              % nearest integer y-coordinate
ix = round(v(1,1));              % nearest integer x-coordinate

iobs = get(hAxes,'Userdata');    % indices of currently displayed obs
i = iobs(iy);                    % this is the index of the selected obs

% create a window with descriptive information

obsinfo(odssubset(ods,i))

%-------------------------------------------------------------------------

function clickonnum

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)

% in case of right-click, do nothing:

if strcmp(get(hFig,'SelectionType'),'alt'), return, end

% get mouse click position and index of the corresponding obs:

v = get(findobj(hFig,'Tag','numaxes'),'CurrentPoint');
iy = round(v(1,2));              % nearest integer y-coordinate

hAxes = findobj(hFig,'Tag','listaxes');  % handle to axes containing obs list
iobs = get(hAxes,'Userdata');    % indices of currently displayed obs
i = iobs(iy);                    % this is the index of the selected obs

% get ods stored with the list:

ods = get(hFig,'Userdata');

% toggle selection flag:

ods.select(i) = ~ods.select(i);

% store modified ods with the list:

set(hFig,'Userdata',ods)

% re-display data

dsplay(iobs(1))

%-------------------------------------------------------------------------
function obsinfo(ods)

[KTS,KXS,QCHS,QCXS,KTPRS,KTSFC,KTRAD] = ...
    dconfig('KTS','KXS','QCHS','QCXS','KTPRS','KTSFC','KTRAD');

if isfield(ods,'cidx'), 
    colr = ods.cinfo(ods.cidx).rgb; 
else
    colr = [0.8125 1.0000 0.1875];
end

isprs = any(ods.kt==KTPRS);
issfc = any(ods.kt==KTSFC);
israd = any(ods.kt==KTRAD);
if issfc, levstr = 'surface';
elseif israd, levstr = int2str(ods.lev);
elseif isprs, levstr = [num2str(ods.lev) ' Pa']; 
elseif ods.kt==22, levstr = [num2str(ods.lev) ' Pa']; 
else, levstr = 'undefined'; end

lonlat = [lonstr(ods.lon) ', ' latstr(ods.lat)];
locstr = [lonlat ', ' levstr];
timstr = jdaystr(ods.first_julian_day + double(ods.time)/1440,0);

units  = KTS(ods.kt==[KTS.value]).units;
obsstr = [num2str(ods.obs) ' ' units];

if isfinite(ods.omf),
    omfstr = [num2str(ods.omf) ' ' units];
    fstr = [num2str(ods.obs-ods.omf) ' ' units];
else,
    omfstr = 'undefined';
    fstr = 'undefined';
end

if isfinite(ods.oma),
    omastr = [num2str(ods.oma) ' ' units];
    astr = [num2str(ods.obs-ods.oma) ' ' units];
else,
    omastr = 'undefined';
    astr = 'undefined';
end

if isfinite(ods.xm), u1str = num2str(ods.xm);
else, u1str = 'undefined'; end

if isfinite(ods.sigo), u2str = [num2str(ods.sigo) ' ' units];
else, u2str = 'undefined'; end

txt{ 1,1} = 'Data source:';              txt{ 1,2} = KXS(ods.kx==[KXS.value]).id;
txt{ 2,1} = 'Sounding index:';           txt{ 2,2} = int2str(ods.ks);
txt{ 3,1} = 'Data type:';                txt{ 3,2} = KTS(ods.kt==[KTS.value]).id;
txt{ 4,1} = 'Observation time:';         txt{ 4,2} = timstr;
txt{ 5,1} = 'Latitude:';                 txt{ 5,2} = latstr(ods.lat);
txt{ 6,1} = 'Longitude:';                txt{ 6,2} = lonstr(ods.lon);
txt{ 7,1} = 'Level:';                    txt{ 7,2} = levstr;
if israd, txt{ 7,1} = 'Channel number:'; end
txt{ 8,1} = 'Observed value:';           txt{ 8,2} = obsstr;
txt{ 9,1} = 'Observation error:';        txt{ 9,2} = u2str;
txt{10,1} = 'Background value:';         txt{10,2} = fstr;
txt{11,1} = 'Background residual:';      txt{11,2} = omfstr;
txt{12,1} = 'Analyzed value:';           txt{12,2} = astr;
txt{13,1} = 'Analysis residual:';        txt{13,2} = omastr;
txt{14,1} = 'QC history flag:';          txt{14,2} = QCHS(ods.qch==[QCHS.value]).id;
txt{15,1} = 'QC rejection flag:';        txt{15,2} = QCXS(ods.qcx==[QCXS.value]).id;
txt{16,1} = 'Metadata:';                 txt{16,2} = u1str;

nlines = size(txt,1);
l1 = 1; for j = 1:nlines, l1 = max(l1,length(txt{j,1})); end
l2 = 1; for j = 1:nlines, l2 = max(l2,length(txt{j,2})); end
i1 = 2; i2 = 30; imax = 60; jmax = nlines + 2;

hFig = figure('Units','character','Position',[60 20 imax jmax], ...
    'MenuBar','none','NumberTitle','off','Resize','off', ...
    'Name',['Observation at ' locstr], ...
    'Visible','on','Tag','obsinfo');

hAx = axes('Parent',hFig,'Units','character','Position',[0 0 imax jmax], ...
    'Box','on','DrawMode','fast','color',[0 0 0], ...
    'Xtick',-1,'Ytick',-1,'XtickLabel','','YtickLabel','', ...
    'XLim',[1 imax],'YLim',[1 jmax],'YDir','reverse');

for j = 1:nlines,
    text(i1,j+1,txt{j,1},'Parent',hAx,'FontName','default','Color',[1 1 1]);
    text(i2,j+1,txt{j,2},'Parent',hAx,'FontName','default','Color',colr);
end

%-------------------------------------------------------------------------
function reorder(option)

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)
ods = get(hFig,'Userdata');

hOrder = findobj(hFig,'Tag','descending');

descending = strcmp(get(hOrder,'Checked'),'on');
if descending, order = 'descend'; else, order = 'ascend'; end

switch option

    case 'selected'
        
        [x,i] = sort(ods.select,order);
        ods = odssubset(ods,i);

    case 'unselected'
        
        [x,i] = sort(~ods.select,order);
        ods = odssubset(ods,i);

    case 'descending'
        
        if descending, set(hOrder,'Checked','off'); else, set(hOrder,'Checked','on'); end

    otherwise

        attr = option;
        [x,i] = sort(ods.(attr),order);
        ods = odssubset(ods,i);

end

set(hFig,'Userdata',ods)

dsplay(1)

%-------------------------------------------------------------------------
function psoundings

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)
ods = get(hFig,'Userdata');
            
name  = 'Plot soundings';
KTS   = dconfig('KTS');
attrs = {'obs','omf','oma','xm','sigo'};
ATTRS = {'observed value (obs)','background residual (omf)',...
        'analysis residual (oma)',...
        'metadata (xm)','observation error (sigo)'};

kt = [];                % find all kt for which there are soundings
for kts = sort(unique(ods.kt))',
    %     for t = unique(ods.time)',
    %         i = find(ods.time==t & ods.kt==kts);
    i = find(ods.kt==kts);
    if any(diff(sort(ods.ks(i)))==0), kt = [kt kts]; end
    %     end
end
kt = unique(kt);

if isempty(kt), warndlg('No univariate soundings, sorry.'), return, end

if length(kt)>1,        % restrict to a single kt (for now)
    prompt = 'Select a single data type:';
    ikts = find(ismember([KTS.value],kt));
    i = listdlg('Name',name,'PromptString',prompt,'SelectionMode','single',...
        'ListSize',[200 200],'ListString',{KTS(ikts).id});
    if isempty(i), return; end
    kt = KTS(ikts(i)).value;
    ods = odssubset(ods,ods.kt==kt);
end

prompt = 'Select a single data attribute:';
ix = listdlg('Name',name,'PromptString',prompt,'SelectionMode','single',...
    'ListSize',[200 200],'ListString',ATTRS);
if isempty(ix), return; end

x = ods.(attrs{ix});
y = ods.lev;

hProfFig = figure;
hProfAx = axes('Parent',hProfFig,'Position',[0.1 0.1 0.3 0.8]);

ia = [];
% ts = sort(unique(ods.time)');
% t1 = min(ts); t2 = max(ts);
% cmap = cool; nc = size(cmap,1);
colr = 'k';
% for t = ts, 
    for ks = unique(ods.ks)';
%         i = find(ods.time==t & ods.ks==ks);
        i = find(ods.ks==ks);
        if length(i)>1,
            h = line(x(i),y(i)); 
%             if diff([t1 t2]),   % multiple times
%                 colr = cmap(floor(1+(nc-1)*(t-t1)/(t2-t1)),:);
%             end
            set(h,'LineStyle','-','LineWidth',1,'Color',colr)
            ia = [ia; i];
        end
    end
% end

ods = odssubset(ods,ia);         % all plotted data
h = line(x(i),y(i)); 
set(h,'LineStyle','none','Marker','.','Markersize',16)

y = unique(ods.lev);
ytick = y; while length(ytick)>16, ytick = ytick(1:2:end); end
set(hProfAx,'PlotBoxAspectRatio',[1 2 1],...
    'YScale','log','YDir','reverse',...
    'YLim',[min(y) max(y)],'YTick',ytick,'YMinorTick','off',...
    'XGrid','on','YGrid','on','YMinorGrid','off')
title({KTS(kt==[KTS.value]).id; ATTRS{ix}})        

%             hLeg = legend(choices{ix});
%             pos = get(hLeg,'Position'); pos(1) = 0.5, set(hLeg,'Position',pos)

% Create context menus for the plot

hMenu = uicontextmenu;
set(gca,'UIContextMenu',hMenu);

uimenu(hMenu,'Label','Linear scale',...
    'Callback','set(gca,''Yscale'',''linear'')');
uimenu(hMenu,'Label','Logarithmic scale',...
    'Callback','set(gca,''Yscale'',''log'')');


%-------------------------------------------------------------------------
function ptseries

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)
ods = get(hFig,'Userdata');

name   = 'Plot timeseries';
prompt = 'Select one or more data attributes:';
iy = listdlg('Name',name,'PromptString',prompt,...
    'ListSize',[200 200],'ListString',choices);
if isempty(iy), return; end
%     x = datenum(jdaystr(ods.first_julian_day + ods.time/1440,0)); 
%     datetick('x',1)


%-------------------------------------------------------------------------
function mkmenubar(hFig)

ATTRS = dconfig('OBSATTRIBUTES');

% File menu 

h = uimenu(hFig,'Label','File'); 
uimenu(h,'Label','Print...',...
    'Callback','printdlg');
uimenu(h,'Label','Save...',...
    'Callback','filemenufcn FileSaveAs');
uimenu(h,'Label','Page Setup...',...
    'Callback','filemenufcn FilePageSetup');
uimenu(h,'Label','Close','Separator','on',...
    'Callback','filemenufcn FileClose');

% Sort menu

h = uimenu(hFig,'Label','Sort');
for i = 1:length(ATTRS.names),
    uimenu(h,'Label',['By ' ATTRS.names{i} ' (' ATTRS.descr{i} ')'],...
        'Callback',['obslist reorder ' ATTRS.names{i}])
end
uimenu(h,'Separator','on',...
    'Label','Selected on top','Callback','obslist reorder selected');
uimenu(h,...
    'Label','Unselected on top','Callback','obslist reorder unselected');
uimenu(h,'Separator','on',...
    'Tag','descending','Label','Sort in descending order','Checked','off',...
    'Callback','obslist reorder descending');

% Select menu

h = uimenu(hFig,'Label','Select'); 
hs = uimenu(h,'Label','Based on attribute values...');
for i = 1:length(ATTRS.names),
    uimenu(hs,'Label',[ATTRS.names{i} ' (' ATTRS.descr{i} ')'],...
        'Callback',['obslist select attr ' ATTRS.names{i}]);
end
uimenu(h,'Separator','on',...
    'Label','Select all','Callback','obslist select all')
uimenu(h,'Label','Unselect all','Callback','obslist select none')
uimenu(h,'Separator','on',...
    'Label','Create a new list with selected data only',...
    'Callback','obslist select list selected')
uimenu(h,'Label','Create a new list with unselected data only',...
    'Callback','obslist select list unselected')
uimenu(h,'Separator','on',...
    'Label','Send selected data to workspace',...
    'Callback','ods=get(gcbf,''Userdata''); ods=odssubset(ods,ods.select)')
uimenu(h,'Label','Send unselected data to workspace',...
    'Callback','ods=get(gcbf,''Userdata''); ods=odssubset(ods,~ods.select)')

% Plot menu

h = uimenu(hFig,'Label','Plot'); 
hs = uimenu(h,'Label','Create an interactive map...');
uimenu(hs,'Label','using selected data only',...
    'Callback','ods=get(gcbf,''Userdata''); obsmap(odssubset(ods,ods.select))')
uimenu(hs,'Label','using unselected data only',...
    'Callback','ods=get(gcbf,''Userdata''); obsmap(odssubset(ods,~ods.select))')
uimenu(hs,'Label','using all data',...
    'Callback','obsmap(get(gcbf,''Userdata''))')   
hs = uimenu(h,'Label','Create a summary plot...');
uimenu(hs,'Label','using selected data only',...
    'Callback','ods=get(gcbf,''Userdata''); obsplot(odssubset(ods,ods.select))')
uimenu(hs,'Label','using unselected data only',...
    'Callback','ods=get(gcbf,''Userdata''); obsplot(odssubset(ods,~ods.select))')
uimenu(hs,'Label','using all data',...
    'Callback','obsplot(get(gcbf,''Userdata''))')   
hs = uimenu(h,'Label','Create a histogram...');
uimenu(hs,'Label','using selected data only',...
    'Callback','ods=get(gcbf,''Userdata''); obshist(odssubset(ods,ods.select))')
uimenu(hs,'Label','using unselected data only',...
    'Callback','ods=get(gcbf,''Userdata''); obshist(odssubset(ods,~ods.select))')
uimenu(hs,'Label','using all data',...
    'Callback','obshist(get(gcbf,''Userdata''))')  

% Options menu

h = uimenu(hFig,'Label','Options');
uimenu(h,'Label','Set background color',...
    'Callback','set(gcf,''Color'',uisetcolor(gcf))')

% Help menu

h = uimenu(hFig,'Label','Help');
uimenu(h,'Label','Basics','Callback','obslist helpme basics') 
uimenu(h,'Separator','on',...
    'Label','Displaying a single observation','Callback','obslist helpme obsinfo') 
if any(findobj(hFig,'Style','slider')), enable = 'on'; else, enable = 'off'; end
uimenu(h,'Label','Scrolling the list','Callback','obslist helpme scroll','Enable',enable) 
uimenu(h,'Label','Selecting observations','Callback','obslist helpme select') 

% enable = 'off'; 
% if nobs<10000,
%     for t = unique(ods.time)',  % see if there are any soundings
%         i = find(ods.time==t);
%         if any(diff(sort(ods.ks(i)))==0), enable = 'on'; end
%     end
% end
% uimenu(hMenu,'Label','Plot soundings..','Enable',enable,...
%     'Callback','obslist psoundings','Separator','on');
% uimenu(hMenu,'Label','Plot timeseries..','Enable','off',...
%     'Callback','obslist ptseries');

%-------------------------------------------------------------------------
function helpme(option)

label = ['Help: ' get(gcbo,'Label')];

switch option
    
    case 'basics'
        
        helpdlg({['This figure displays a complete list of data attributes. Each row '...
            'corresponds to a single observation. If there are more observations '...
            'than fit on a single page, a scrollbar appears on the right edge of the window. '...
            'Observations can be rearranged using the functions in the ''Sort'' menu. '...
            'There are many ways to select groups of observations, which '...
            'can then be used to create a new list, and/or various '...
            'kinds of plots.']},label)
    
    case 'obsinfo'
        
        helpdlg({['Click on any number in the list to '...
           'display descriptive labels for all the attributes']},label)
       
    case 'select'
        
        helpdlg({'Select (or unselect) observations using the following methods: ';...
           ' ';...
           '- click on the sequence numbers in the leftmost column of the list;';...
           '- click on the buttons at the top of the columns; or';... 
           '- use one of the functions in the ''Select'' menu.'},label)
       
    case 'scroll'
        
        helpdlg({'Scroll using the slider control on the right edge of the window: ';...
           ' ';...
           '- drag the slider bar up or down';...
           '- click on the arrows at top or bottom; or';... 
           '- click anywhere else on the slider.'},label)
       
end

