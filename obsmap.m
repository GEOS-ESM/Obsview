function obsmap(arg,varargin)

% OBSMAP Plot interactive observation map.
%
% OBSMAP(ODS) produces an interactive map showing the observations in the
%     ods structure ODS.

% 19Oct2001 Dick Dee (dee@dao.gsfc.nasa.gov)
% 23Apr2004 Dick Dee - ECMWF version

% check whether this is a callback or not

if ischar(arg), % must be a callback
    feval(arg,varargin{:})      % simply pass to the callback routine..
    return                      % .. and get out.
else            % otherwise this will set up a new plot
    ods = arg;  % first arg is always the ods structure with the data
end

% check if we have data to plot

ndata = length(ods.kt);
if ndata==0, return, end

% confirm for many data

if ndata>5000,
    go = questdlg({['Your data selection contains ' int2str(ndata) ' observations.'] ...
            ' ' ...
            'Plotting this map may take a while. Are you sure you want' ...
            'to continue?'},'Just checking..','Yes','No','No');
    if strcmp(go,'No'), return, end
end

% get handle to obsview gui, if it exists

hGui = findobj('Tag','ObsviewGui');  

% decide whether to use map projections

if isempty(ver('map')),    % need the Mapping Toolbox
    mapproj = 0;  
else, 
    mapproj = 1; 
    if ~isempty(hGui) && strcmp(get(findobj(hGui,'Tag','mapproj'),'Checked'),'off'),
        mapproj = 0;                   % turned off by obsview gui 
    end
end

% get central longitude for global maps

centlon = 0;
if ~isempty(hGui) && strcmp(get(findobj(hGui,'Tag','dateline'),'Checked'),'on'), 
    centlon = 180; 
end

% create a new figure window

hFig = figure;
set(hFig,'Tag','obsmap','Units','pixels')

% store data and options with the map figure

set(hFig,'Userdata',ods)
setappdata(hFig,'mapproj',mapproj)
setappdata(hFig,'centlon',centlon)

% create the different components of the figure

if mapproj, dx = 0.02; dy = 0.02; 
else,       dx = 0.05; dy = 0.05; end
dt = 0.10; 
MapPosition = [dx dy 1-2*dx 1-4*dy-dt];
TtlPosition = [dx 1-dy-dt 1-2*dx dt];

% create a title

figure(hFig)
axes('Parent',hFig,'Tag','title','Position',TtlPosition,'NextPlot','add',...
    'Visible','off');
ptitle('create')

% create the map

axes('Parent',hFig,'Tag','map','Position',MapPosition,'NextPlot','add');
pmap('create')


%-------------------------------------------------------------------------
function pmap(action,varargin)

switch action

case 'create'

    hAxes = gca;
    hFig = gcf;
    ods = get(hFig,'Userdata');
    mapproj = getappdata(hFig,'mapproj');  % see if using a map projection
    centlon = getappdata(hFig,'centlon');  % central longitude for global maps
    
    if isfield(ods,'sslon')&&any(ods.sslon), xlim = ods.sslon; else, xlim = [-180 180]; end
    if isfield(ods,'sslat')&&any(ods.sslat), ylim = ods.sslat; else, ylim = [ -90  90]; end
    
    % determine longitude limits xlim such that:
    %                -360<xlim(1)<xlim(2)<360 and xlim(2)-xlim(1)<=360
    
    xlim = -180 + mod(xlim + 180,360);     % xlim in [-180,180)
    
    if xlim(1)==xlim(2),                   % this means all longitudes: global map
        xlim = centlon + [-180 180];
    elseif 0<=xlim(2) && xlim(2)<xlim(1),
        xlim(1) = xlim(1) - 360;
    elseif xlim(2)<0  && xlim(2)<xlim(1),
        xlim(2) = xlim(2) + 360;
    end
    
    if mapproj,    % using map projection 
        
        setupprojection(ylim,xlim)   % creates a suitable map projection depending on the lon-lat limits
        set(gca,'Visible','off')
        setm(gca,'Grid','on',...
            'meridianlabel','on','parallellabel','on',...
            'GLineStyle','--','LabelFormat','none',...
            'frame','on','FLineWidth',1,'FFaceColor',[1 1 1])
        
        load coast
        plotm(lat,long,'Tag','coast','Color',0.3*[1 1 1])
        
        [x,y] = mfwdtran(double(ods.lat),double(ods.lon)); 
        
    else,          % using rectangular lon-lat map
        
        set(hAxes,'XLim',xlim,'YLim',ylim,'Box','on','XGrid','on','YGrid','on');
        if diff(xlim)==360, set(hAxes,'XTick',linspace(xlim(1),xlim(2),7)); end
        if diff(ylim)==180, set(hAxes,'YTick',linspace(ylim(1),ylim(2),7)); end
        
        load coast
        plot(long,lat,'Color',0.3*[1 1 1])   % brute-force dateline crossing
        if xlim(1)<-180, plot(long-360,lat,'Color',0.3*[1 1 1]); end 
        if xlim(2)> 180, plot(long+360,lat,'Color',0.3*[1 1 1]); end 
        
        % convert ods longitudes to map lonitudes, if necessary
        
        x = ods.lon; 
        y = ods.lat;
        if xlim(1)<-180, i = (x>xlim(1)+360); x(i) = x(i) - 360; end
        if xlim(2)> 180, i = (x<xlim(2)-360); x(i) = x(i) + 360; end
        
    end
    
    % draw the patches
    
    ods.hPatch = NaN*ones(size(ods.lon));
    for i = 1:length(x),
        ods.hPatch(i) = patch('Parent',hAxes,'Vertices',[x(i) y(i) 1],'Faces',1, ...
            'FaceVertexCData',NaN,'FaceColor','none','EdgeColor','none', ...
            'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor','k','Tag','obspatch');
    end
    
    % store ods with patch handles with the figure
    
    set(hFig,'Userdata',ods)

    % set the title

    title('Observation Locations','Color','blue')
    
    % check if any windvectors
    
    KTWND = dconfig('KTWND');
    iswnd = 0; 
    for i = 1:length(KTWND),
        ktwnd = KTWND{i};
        iswnd = iswnd | (any(ods.kt==ktwnd(1)) & any(ods.kt==ktwnd(2)));
    end

    % create the contextmenu, and store its handle with the map figure
    
    hMenu = uicontextmenu;
    setappdata(hFig,'contextmenu',hMenu)
    
    % attach to the axes and its children
    
    set([hAxes; allchild(hAxes)],'UIContextMenu',hMenu);
        
    hColor = uimenu(hMenu,'Label','Color locations by');
    uimenu(hColor,'Label','data type index (kt)',                         'Callback','obsmap pmap color kt')
    uimenu(hColor,'Label','data source index (kx)',                       'Callback','obsmap pmap color kx')
    if any(isfinite(ods.lev)),        enable = 'on'; else, enable = 'off'; end
    uimenu(hColor,'Label','level',                        'Enable',enable,'Callback','obsmap pmap color lev')
    uimenu(hColor,'Label','observation time (time)',                      'Callback','obsmap pmap color time')
    uimenu(hColor,'Label','observed value (obs)',                         'Callback','obsmap pmap color obs')
    if any(isfinite(ods.sigo)),       enable = 'on'; else, enable = 'off'; end
    uimenu(hColor,'Label','observation error (sigo)',     'Enable',enable,'Callback','obsmap pmap color sigo')
    if any(isfinite(ods.omf)),        enable = 'on'; else, enable = 'off'; end
    uimenu(hColor,'Label','background residual (omf)',    'Enable',enable,'Callback','obsmap pmap color omf')
    if any(isfinite(ods.oma)),        enable = 'on'; else, enable = 'off'; end
    uimenu(hColor,'Label','analysis residual (oma)',      'Enable',enable,'Callback','obsmap pmap color oma')
    uimenu(hColor,'Label','QC exclusion mark (qcx)',                      'Callback','obsmap pmap color qcx')
    uimenu(hColor,'Label','QC history mark (qch)',                        'Callback','obsmap pmap color qch')
    if any(isfinite(ods.xm)),         enable = 'on'; else, enable = 'off'; end
    uimenu(hColor,'Label','metadata (xm)',                'Enable',enable,'Callback','obsmap pmap color xm')
    
    if iswnd, enable = 'on'; else, enable = 'off'; end
    uimenu(hColor,'Label','observed wind speed (ws_obs)','Separator','on','Enable',enable, ...
        'Callback','obsmap pmap color ws_obs')
    if ~any(isfinite(ods.omf)), enable = 'off'; end
    uimenu(hColor,'Label','wind speed background residual (ws_omf)','Enable',enable, ...
        'Callback','obsmap pmap color ws_omf')
    if ~any(isfinite(ods.oma)), enable = 'off'; end
    uimenu(hColor,'Label','wind speed analysis residual (ws_oma)','Enable',enable, ...
        'Callback','obsmap pmap color ws_oma')
    
    hCmap = uimenu(hMenu,'Label','Color scheme','Tag','colorscheme','Enable','off');
    uimenu(hCmap,'Label','Blue-white-red',...
        'Callback','obsmap pmap cmap bwr')
    uimenu(hCmap,'Label','Blue-cyan-yellow-orange-red',...
        'Callback','obsmap pmap cmap jet')
    uimenu(hCmap,'Label','Cyan-magenta',...
        'Callback','obsmap pmap cmap cool')
    
    uimenu(hMenu,'Label','Set color scale','Tag','cscale','Enable','off',...
        'Callback','obsmap plegend cscale')
    
    if iswnd, enable = 'on'; else, enable = 'off'; end
    hWind = uimenu(hMenu,'Label','Show wind vectors','Separator','on','Enable',enable);
    uimenu(hWind,'Label','observed','Callback','obsmap pwind create obs')
    if any(isfinite(ods.omf)), enable = 'on'; else, enable = 'off'; end
    uimenu(hWind,'Label','background residual','Enable',enable,...
        'Callback','obsmap pwind create omf')
    if any(isfinite(ods.oma)), enable = 'on'; else, enable = 'off'; end
    uimenu(hWind,'Label','analysis residual','Enable',enable,...
        'Callback','obsmap pwind create oma')
    if any(findobj(hAxes,'Tag','pwind')), enable = 'on'; else, enable = 'off'; end
    uimenu(hMenu,'Label','Remove wind vectors','Tag','delwind','Enable',enable,...
        'Callback','obsmap pwind delete')
    
    % uimenu(hMenu,'Label','Add scalar field','Separator','on',...
    %     'Callback','obsmap pfield create')
    
    uimenu(hMenu,'Label','Create a summary plot',...
        'Callback','obsplot(get(gcbf,''Userdata''))','Separator','on');
    uimenu(hMenu,'Label','Create a list of attribute values',...
        'Callback','obslist(get(gcbf,''Userdata''))')
    uimenu(hMenu,'Label','Create a histogram','Tag','histogram',...
        'Callback','obshist(get(gcbf,''Userdata''),getappdata(gcbf,''coloring''))')
    uimenu(hMenu,'Label','Put ods in workspace','Separator','on',...
        'Callback','ods=get(gcbf,''Userdata'')');
    uimenu(hMenu,'Label','Duplicate this figure',...
        'Callback','copyobj(gcbf,0)')
    
    mapproj = getappdata(hFig,'mapproj');
    if mapproj, enable = 'on'; else, enable = 'off'; end
    hMap = uimenu(hMenu,'Label','Map appearance','Enable',enable,'Separator','on');
    if mapproj, 
        checked = get(findobj(hAxes,'Tag','coast'),'Visible');    
        uimenu(hMap,'Label','Show coastline','Tag','coastline','Checked',checked,...
            'Callback','obsmap pmap toggle coastline')
        checked = getm(hAxes,'Grid');
        uimenu(hMap,'Label','Show map grid','Tag','mapgrid','Checked',checked,...
            'Callback','obsmap pmap toggle mapgrid')
        checked = getm(hAxes,'MeridianLabel');
        uimenu(hMap,'Label','Label meridians','Tag','mlabels','Checked',checked,...
            'Callback','obsmap pmap toggle mlabels')
        checked = getm(hAxes,'ParallelLabel');
        uimenu(hMap,'Label','Label parallels','Tag','plabels','Checked',checked,...
            'Callback','obsmap pmap toggle plabels')
    end
    
    % prepare callback

    set([hAxes; allchild(hAxes)],'ButtonDownFcn','obsmap pobslist')

case 'color'

    attr = varargin{1};

    [h,hFig] = gcbo;                     % handles to cbo (h) and figure (hFig)
    hAxes = findobj(hFig,'Type','axes','Tag','map');   % handle to map axes

    watchon;                             % this might take a while

    ods = get(hFig,'Userdata');          % get ods stored with the map

    if ~isfield(ods,attr),               % new attributes to be defined
        switch attr
        case {'ws_obs','ws_omf','ws_oma'},  % related to wind speed
            
            % indices of matching wind components
            
            if ~(isfield(ods,'iu')&&isfield(ods,'iv')),
                [ods.iu,ods.iv] = uvmatch(ods); 
            end
            
            ods.ws_obs = NaN*ones(size(ods.obs));
            ods.ws_omf = NaN*ones(size(ods.obs));
            ods.ws_oma = NaN*ones(size(ods.obs));
            
            ws_obs = sqrt( ods.obs(ods.iu)                 .^2 ...
                        +  ods.obs(ods.iv)                 .^2);
            ws_fcs = sqrt((ods.obs(ods.iu)-ods.omf(ods.iu)).^2 ...
                        + (ods.obs(ods.iv)-ods.omf(ods.iv)).^2);
            ws_ana = sqrt((ods.obs(ods.iu)-ods.oma(ods.iu)).^2 ...
                        + (ods.obs(ods.iv)-ods.oma(ods.iv)).^2);
            ods.ws_obs(ods.iu) = ws_obs;
            ods.ws_omf(ods.iu) = ws_obs - ws_fcs;
            ods.ws_oma(ods.iu) = ws_obs - ws_ana;
            ods.ws_obs(ods.iv) = ws_obs;
            ods.ws_omf(ods.iv) = ws_obs - ws_fcs;
            ods.ws_oma(ods.iv) = ws_obs - ws_ana;
            
            set(hFig,'Userdata',ods)     % store ods with new attributes
            
        end
    end

    c = ods.(attr);              % now define the color indices

    for i = 1:length(c), set(ods.hPatch(i),'FaceVertexCData',c(i)), end   % set the patch colors

    ptitle('replace',['Colored by ' attr])
    setappdata(hFig,'coloring',attr)                 % store coloring attr
    axes(hAxes)
    plegend('create')                                % create a legend
    set(findobj('Tag','colorscheme'),'Enable','on')  % allow user to change the color scheme

    watchoff;                            % done

case 'cmap'

    scheme = varargin{1};

    switch scheme
    case 'bwr'         % create a blue-white-red colormap
        a = linspace(0,1,32)';
        b = linspace(1,0,32)';
        c = ones(32,1);
        colormap([a a c; c b b])
    otherwise
        colormap(scheme)
    end
    
case 'toggle'
    
    hFig = gcf; 
    hAxes = findobj(hFig,'Type','axes','Tag','map');   % handle to map axes
    
    option = varargin{1};
    hMenu = findobj(hFig,'Type','uimenu','Tag',option);
    onoroff = get(hMenu,'Checked');
    if strcmp(onoroff,'on'), onoroff = 'off'; else, onoroff = 'on'; end
    
    calledsetm = 0;
    switch option
    case 'mapgrid', 
        setm(hAxes,'Grid',onoroff,...
            'MeridianLabel',onoroff,...
            'ParallelLabel',onoroff); 
        set(findobj(hFig,'Type','uimenu','Tag','mapgrid'),'Checked',onoroff)
        set(findobj(hFig,'Type','uimenu','Tag','mlabels'),'Checked',onoroff)
        set(findobj(hFig,'Type','uimenu','Tag','plabels'),'Checked',onoroff)
        calledsetm = 1; 
    case 'mlabels', 
        setm(hAxes,'MeridianLabel',onoroff); 
        set(findobj(hFig,'Type','uimenu','Tag','mlabels'),'Checked',onoroff)
        calledsetm = 1;
    case 'plabels', 
        setm(hAxes,'ParallelLabel',onoroff); 
        set(findobj(hFig,'Type','uimenu','Tag','plabels'),'Checked',onoroff)
        calledsetm = 1; 
    case 'coastline', 
        hCoast = findobj(hAxes,'Tag','coast'); 
        set(hCoast,'Visible',onoroff)
        set(findobj(hFig,'Type','uimenu','Tag','coastline'),'Checked',onoroff)
    end
    
    if calledsetm,  % must reattach contextmenu and 'ButtonDownFcn' callback
        hCMenu = getappdata(hFig,'contextmenu');
        set([hAxes; allchild(hAxes)],'UIContextMenu',hCMenu);
        set([hAxes; allchild(hAxes)],'ButtonDownFcn','obsmap pobslist create')
    end
    
end

%-------------------------------------------------------------------------
function plegend(action,varargin)

[h,hFig] = gcbo;                         % handles to cbo (h) and figure (hFig)
hAxes = findobj(hFig,'Type','axes','Tag','map');       % handle to map axes

attr = getappdata(hFig,'coloring');      % get name of coloring attribute
ods = get(hFig,'Userdata');
c = double(ods.(attr));          % get color indices

switch action

case 'create'

    if nargin==2, cscale = varargin{1}; else, cscale = []; end

    cf = c(~isnan(c));        % finite color indices
    cs = unique(cf);          % unique finite color indices..
    cs = -sort(-cs(:));       % ..column vector in descending order

    switch attr
    case {'kx','kt','qcx','qch'},
        continuous = 0;                   % discrete range of colors
        cscale = cs;                      % use all colors for legend
    otherwise
        continuous = 1;                   % continuous range of colors
        if isempty(cscale),               % create colorscale for legend..
            if ~isempty(cf),
                if length(cs)>10,         % ..use 7 linearly spaced colors
                    cscale = linspace(max(cf),min(cf),7)';
                else
                    cscale = cs;          % ..use all colors
                end
            else
                cscale = [];              % ..only NaNs
            end
        end
    end

    legend off                            % get rid of the old legend, if any

    % set the color scale for this plot

    if length(cscale)>1,       set(hAxes,'CLim',[min(cscale) max(cscale)]);
    elseif length(cscale)==1,  set(hAxes,'CLim',[cscale-1    cscale+1   ]);  % no range
    else              % no data
    end

    % allow user to modify the color scale for continuous attributes

    if continuous, enable = 'on'; else, enable = 'off'; end
    set(findobj('Tag','cscale'),'Enable',enable)

    % create the legend

    if any(isnan(c)), cscale = [NaN; cscale]; end
    for i = 1:length(cscale),
        ch(i) = patch([0 1 1 0],[1 1 0 0],cscale(i),'Visible','off');  % create a patch for each color
    end
    legend(ch,num2str(cscale),3);     % draw the new legend

case 'cscale'

    prompt = ['Enter new legend values: ' ...
            '(data range from ' num2str(min(c)) ' to ' num2str(max(c)) ')'];
    answer = inputdlg(prompt,'Set color scale',1);
    if isempty(answer), return, end
    cscale = eval(['[' answer{1} ']']);
    if length(cscale)>30, errordlg('Too many values; try again.'); end
    cscale = -sort(-cscale(:));  % column vector in descending order

    plegend('create',cscale)   % recreate the legend

end

%-------------------------------------------------------------------------
function pwind(action,varargin)

[h,hFig] = gcbo;         % handles to cbo (h) and figure (hFig)
hAxes = findobj(hFig,'Type','axes','Tag','map');      % handle to map axes

switch action

case 'create'
        
    attr = varargin{1};
        
    ods = get(hFig,'Userdata');
    
    % find indices of all matching wind components
    
    if ~(isfield(ods,'iu')&&isfield(ods,'iv')),
        [ods.iu,ods.iv] = uvmatch(ods); 
        set(hFig,'Userdata',ods)
    end
       
    switch attr
    case 'obs', 
        u = double(ods.obs(ods.iu));
        v = double(ods.obs(ods.iv));
    case 'omf', 
        u = double(ods.omf(ods.iu));
        v = double(ods.omf(ods.iv));
    case 'oma', 
        u = double(ods.oma(ods.iu));
        v = double(ods.oma(ods.iv));
    end
    lon = double(ods.lon(ods.iu));
    lat = double(ods.lat(ods.iu));
    
    delete(findobj(hAxes,'Tag','pwind'))   % get rid of the old wind vectors, if any
    
    if getappdata(hFig,'mapproj'),   % using a map projection
        
        dlat = v;
        dlon = u;
        i = (abs(lat)<90); dlon(i) = u(i)./cos(pi*lat(i)/180);
        hWind = quiverm(lat,lon,dlat,dlon,2);    % draw the arrows
        
    else                             % using a rectangular map
        
        xlim = get(hAxes,'XLim');    % convert ods longitudes to map longitudes, if necessary
        if xlim(1)<-180, i = (lon>xlim(1)+360); lon(i) = lon(i) - 360; end
        if xlim(2)> 180, i = (lon<xlim(2)-360); lon(i) = lon(i) + 360; end
        
        hWind = quiver(lon,lat,u,v);     % draw the arrows
        
    end
    
    set(hWind,'Tag','pwind','Hittest','off','LineWidth',1,'color','b');
    set(findobj(hFig,'Tag','delwind'),'Enable','on')  % allow user to delete
    
case 'delete'

    delete(findobj(hAxes,'Tag','pwind'))
    set(findobj(hFig,'Tag','delwind'),'Enable','off')

end

%-------------------------------------------------------------------------
function pobslist

[h,hFig] = gcbo;                    % handles to cbo (h) and figure (hFig)
hAxes = findobj(hFig,'Type','axes','Tag','map');  % handle to map axes

% in case zoom is on, do nothing:

if isstruct(getappdata(hFig,'ZOOMFigureState')), return, end

% in case of right-click, do nothing:

if strcmp(get(hFig,'SelectionType'),'alt'), return, end

% get ods stored with the map:

ods = get(hFig,'Userdata');

% get mouse click and release positions:

va = get(hAxes,'CurrentPoint');
x0 = va(1,1);
y0 = va(1,2);
rbbox;
va = get(hAxes,'CurrentPoint');
x1 = va(1,1);
y1 = va(1,2);

mapproj = getappdata(hFig,'mapproj'); 
if mapproj,               % using a map projection
    
    [lat0,lon0] = minvtran(x0,y0);
    [lat1,lon1] = minvtran(x1,y1);
    lon = ods.lon;        % convert ods longitudes to map longitudes, if necessary
    
else,                     % using a rectangular map
    
    lat0 = y0; lon0 = x0;
    lat1 = y1; lon1 = x1;
    lon = ods.lon;        % convert ods longitudes to map longitudes, if necessary
    xlim = get(hAxes,'XLim');
    if xlim(1)<-180, i = (lon>xlim(1)+360); lon(i) = lon(i) - 360; end
    if xlim(2)> 180, i = (lon<xlim(2)-360); lon(i) = lon(i) + 360; end
    
end

if (lat1==lat0 && lon1==lon0),   % find data at location closest to the mouse click:
    
    c = cos(pi*lat0/180);
    dmin =  min(c*abs(lon0-lon)+abs(lat0-ods.lat));
    if dmin>3, return; end
    iobs = find(c*abs(lon0-lon)+abs(lat0-ods.lat)==dmin);
    
else  % find all data within the rectangle defined by the click/release positions:
    
    iobs = find(min(lon0,lon1)<=    lon &     lon<=max(lon0,lon1) &...
        min(lat0,lat1)<=ods.lat & ods.lat<=max(lat0,lat1));
    
end

nobs = length(iobs);  % number of selected obs

if nobs==0, return; end

[x,i] = sort(ods.lev(iobs)); % sort by level
ods = odssubset(ods,iobs(i));

obslist(ods)

%-------------------------------------------------------------------------
function ptitle(action,varargin)

switch action

case 'create'

    hAxes = gca;
    ods = get(gcf,'Userdata');

    ssid  = []; if isfield(ods,'ssid' ), ssid  = ods.ssid ; end
    sslat = []; if isfield(ods,'sslat'), sslat = ods.sslat; end
    sslon = []; if isfield(ods,'sslon'), sslon = ods.sslon; end
    sslev = []; if isfield(ods,'sslev'), sslev = ods.sslev; end
    sskt  = []; if isfield(ods,'sskt' ), sskt  = ods.sskt ; end
    sskx  = []; if isfield(ods,'sskx' ), sskx  = ods.sskx ; end
    ssqcx = []; if isfield(ods,'ssqcx'), ssqcx = ods.ssqcx; end
    ssqch = []; if isfield(ods,'ssqch'), ssqch = ods.ssqch; end

    nobs = length(ods.kt);                % number of obs

    line1 = [ssid ': ' int2str(nobs) ' observations'];
    if ~isempty(sslat),
        str = sprintf('lat=[%6.1f,%6.1f)',min(sslat),max(sslat));
        str = strrep(str,' ',''); str = strrep(str,'.0','');
    else
        str = 'all lat';
    end
    line2 = [str ';'];
    if ~isempty(sslon),
        str = sprintf('lon=[%6.1f,%6.1f)',min(sslon),max(sslon));
        str = strrep(str,' ',''); str = strrep(str,'.0','');
    else
        str = 'all lon';
    end
    line2 = [line2 ' ' str ';'];
    if ~isempty(sslev),
        if length(sslev)==1||sslev(1)==sslev(2), str = sprintf('lev=%6.1f',sslev(1));
        else, str = sprintf('lev=[%6.1f,%6.1f)',min(sslev),max(sslev)); end
        str = strrep(str,' ',''); str = strrep(str,'.0','');
    else
        str = 'all lev';
    end
    line2 = [line2 ' ' str ';'];
    if ~isempty(sskt),
        if length(sskt)<6, str = ['kt=' sprintf('%d,',sskt)]; str(end) = '';
        else, str = sprintf('kt=%d,..,%d',sskt(1),sskt(end));
        end
    else
        str = 'all kt';
    end
    line2 = [line2 ' ' str ';'];
    if ~isempty(sskx),
        if length(sskx)<12, str = ['kx=' sprintf('%d,',sskx)]; str(end) = '';
        else, str = sprintf('kx=%d,..,%d',sskx(1),sskx(end));
        end
    else
        str = 'all kx';
    end
    line2 = [line2 ' ' str ';'];
    if ~isempty(ssqcx),
        if length(ssqcx)<12, str = ['qcx=' sprintf('%d,',ssqcx)]; str(end) = '';
        else, str = sprintf('qcx=%d,..,%d',ssqcx(1),ssqcx(end));
        end
    else
        str = 'all qcx';
    end
    line2 = [line2 ' ' str ';'];
    if ~isempty(ssqch),
        if length(ssqch)<12, str = ['qch=' sprintf('%d,',ssqch)]; str(end) = '';
        else, str = sprintf('qch=%d,..,%d',ssqch(1),ssqch(end));
        end
    else
        str = 'all qch';
    end
    line2 = [line2 ' ' str];
    line3 = []; if isfield(ods,'filename'), line3 = ods.filename; end

    hFig = get(hAxes,'Parent');
    figure(hFig)
    text(0.5,0.9,line1,'Tag','title1',...
        'Interpreter','none','Color','blue','HorizontalAlignment','center','FontSize',11)
    text(0.5,0.5,line2,'Tag','title2',...
        'Interpreter','none','HorizontalAlignment','center','FontSize',9)
    if ~isempty(line3),
        text(0.5,0.2,line3,'Tag','title3',...
            'Interpreter','none','HorizontalAlignment','center','FontSize',9);
    end
    
case 'replace'
    
    line3 = varargin{1};
    hAxes = findobj(gcf,'Type','axes','Tag','title');
    hText = findobj(hAxes,'Tag','title3');
    if any(hText), delete(hText); end
    axes(hAxes);
    text(0.5,0.1,line3,'Tag','title3',...
        'Interpreter','none','Color','red','HorizontalAlignment','center','FontSize',11);

end

%-------------------------------------------------------------------------
function pfield(action,varargin)

if ~strcmp(action,'create'),
    hGui = findobj('Tag','pfieldGui');
    hFig = getappdata(hGui,'mapfig');
    hAxes = findobj(hFig,'Type','axes','Tag','map');
    gfio = getappdata(hFig,'gfio');
end

switch action

case 'create'

    [h,hFig] = gcbo;                    % handles to cbo (h) and figure (hFig)

    % create the gui window:
    % ---------------------
    hGui = figure('Tag','pfieldGui','IntegerHandle','off', ...
        'Units','characters','Position',[100 5 92 7.5], ...
        'Resize','off','Color',[0.95 0.95 0.8], ...
        'NumberTitle','off','Name','Add scalar field','MenuBar','none');

    % keep handle to map figure with the gui:
    % --------------------------------------
    setappdata(hGui,'mapfig',hFig);

    % set up the menubar:
    % ------------------
    hFile = uimenu(hGui,'Label','File');
    uimenu(hFile,'Label','Open gfio file',...
        'Callback','obsmap pfield opengfio');
    uimenu(hFile,'Label','Exit','Separator','on',...
        'Callback','close(gcf)');

    % hOptions = uimenu(hGui,'Label','Options');
    % uimenu(hOptions,'Tag','pcontour','Label','Plot contours',...
    %     'Callback','obsmap pfield toggle pcontour','Checked','on');
    % uimenu(hOptions,'Tag','pscolor','Label','Plot color image',...
    %     'Callback','obsmap pfield toggle pscolor','Checked','off');

    % set up the controls:
    % -------------------
    uicontrol('Tag','Flist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[2 4 88 2], ...
        'Callback','obsmap pfield getginfo');
    uicontrol('Tag','Tlist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[2 1 26 2]);
    uicontrol('Tag','Vlist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[30 1 30 2], ...
        'Callback','obsmap pfield getlevels');
    uicontrol('Tag','Plist','Style','popupmenu','String',{' '}, ...
        'Units','characters','Position',[62 1 18 2]);
    uicontrol('Style','pushbutton','BackgroundColor',[0.8 0.95 0.95], ...
        'Units','characters','Position',[82 1.4 7.8 1.7], ...
        'String','Plot', ...
        'Callback','obsmap pfield plot; close(gcbf)');

    % see if we can find previous file info:
    % -------------------------------------
    gfio = getappdata(hFig,'gfio');  % in the callback figure..
    if isempty(gfio),
        for h = findobj('Tag','obsmap')'
            gfio = getappdata(h,'gfio');  % or else in any other obsmap figure
            if ~isempty(gfio), break; end
        end
    end

    if ~isempty(gfio),  % use previous file name to initialize the gui controls
        pfield('opengfio',gfio.filename)
    else                % or ask the user for new file name
        pfield('opengfio')
    end

case 'opengfio'

    if nargin==1,       % get gfio file name
        [fname, path] = uigetfile('*.*', 'Select a GFIO file:');
        if fname==0, return, end
        gfile = [path fname];
    else
        gfile = varargin{1};
    end
    if ~isgfiofile(gfile),
        warndlg([gfile ' is not a readable gfio file.'])
        return
    end

    % update gui list of gfio files:
    % ----------------------------
    hFlist = findobj(hGui,'Tag','Flist');
    flist = get(hFlist,'string');
    if isempty(deblank(flist{1})), n = 1; else, n = length(flist) + 1; end
    flist{n} = gfile;
    set(hFlist,'String',flist,'Value',n)

    % get file info and update the rest of the gui for this file:
    % ----------------------------------------------------------

    pfield('getginfo',gfile)

case 'getginfo'

    if nargin==1,   % user changed file selection; get its name first
        hFlist = findobj(hGui,'Tag','Flist');
        flist = get(hFlist,'String');
        gfile = flist{get(hFlist,'Value')};  % file name
    else
        gfile = varargin{1};
    end

    % get header info from file and store with map figure:
    % ---------------------------------------------------
    gfio = gfioload(gfile,'info');
    setappdata(hFig,'gfio',gfio)

    hTlist = findobj(hGui,'Tag','Tlist');
    set(hTlist,'String',datestr(gfio.time,0),'Value',1)

    hVlist = findobj(hGui,'Tag','Vlist');
    set(hVlist,'String',[gfio.descr_2D; gfio.descr_3D],'Value',1)

    % set levels:
    % ----------

    pfield('getlevels')

case 'getlevels'

    ivar = get(findobj(hGui,'Tag','Vlist'),'Value');
    hPlist = findobj(hGui,'Tag','Plist');
    if ivar>length(gfio.descr_2D)&&isfield(gfio,'lev'),
        for i = 1:length(gfio.lev),
            str{i} = sprintf('%6.2f%s',gfio.lev(i),gfio.lev_units);
        end
        set(hPlist,'String',str,'Visible','on')
    else
        set(hPlist,'Visible','off')
    end

case 'plot'

    % get file parameters from gui controls

    itime = get(findobj(hGui,'Tag','Tlist'),'Value');     % time index
    ivar = get(findobj(hGui,'Tag','Vlist'),'Value');      % variable index
    if ivar>length(gfio.descr_2D)&&isfield(gfio,'lev'),
        ilev = get(findobj(hGui,'Tag','Plist'),'Value');  % level index, if any
    else
        ilev = [];
    end

    % get the data from file

    gfile = gfio.filename;
    gtime = gfio.time(itime);
    vars  = [gfio.vars_2D; gfio.vars_3D];
    gvar  = vars{ivar};
    if any(ilev),
        glev = gfio.lev(ilev);
        gfio = gfioload(gfile,gvar,'time',gtime,'lev',glev);
    else
        gfio = gfioload(gfile,gvar,'time',gtime);
    end
    z = squeeze(gfio.(gvar)');

    glons = gfio.lon;
    glats = gfio.lat;

    % find grid locations within the map limits

    xlim = get(hAxes,'XLim');
    ylim = get(hAxes,'YLim');

    % convert longitudes to map longitudes, if necessary

    glons = -180 + mod(glons + 180,360);     % glons in [-180,180)
    if xlim(1)<-180, i = (glons>xlim(1)+360); glons(i) = glons(i) - 360; end
    if xlim(2)> 180, i = (glons<xlim(2)-360); glons(i) = glons(i) + 360; end

    i = (ylim(1)-2<=glats & glats<=ylim(2)+2);                 % +/- 2deg extra
    [g,j] = sort(glons); j = j(xlim(1)-2<=g & g<= xlim(2)+2);  % +/- 2deg extra
    glats = glats(i);
    glons = glons(j);
    z     = z(i,j);

    % get new plot parameters from gui controls

%     if strcmp(get(findobj(hGui,'Tag','pcontour'),'Checked'),'on'), pstyle = 'pcontour'; end
%     if strcmp(get(findobj(hGui,'Tag','pscolor'),'Checked'),'on'), pstyle = 'pscolor'; end

    % let's do it

    pstyle = 'pcontour';
    axes(hAxes)

    switch pstyle

    case 'pscolor'  % checkerboard plot

        [lon,lat] = meshgrid(glons,glats);		

        h = pcolor(lon,lat,z);  % creates a single surface object
        set(h,'EdgeColor','interp','FaceColor','interp')

    case 'contourf'   % filled contours

        [lon,lat] = meshgrid(glons,glats);	

        contourf(lon,lat,z);  % creates a bunch of patch objects

    case 'pcontour'    % contours

        [lon,lat] = meshgrid(glons,glats);	

        contour(lon,lat,z,'k');  % creates a bunch of line objects

    end

end

%--------------------------------------------------------------------------
