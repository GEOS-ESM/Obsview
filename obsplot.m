function obsplot(arg,varargin)

% OBSPLOT Plot observation attributes.
%
% OBSPLOT(ODS) produces a plot summarizing the data in the
%     ods structure ODS, and provides various methods for
%     subsetting and replotting the data.

% 25Oct2001 Dick Dee (dee@gmao.gsfc.nasa.gov)
% 23Apr2004 Dick Dee - ECMWF version
% 18Nov2004 Dick Dee - GSI version


if ischar(arg), % must be a callback
    feval(arg,varargin{:})      % simply pass to the callback routine..
    return                      % .. and get out.
else            % otherwise this will set up a new plot
    ods = arg;  % first arg is always the ods structure with the data
    if nargin>1, 
        interactive = false; 
        caller = varargin{1};
    else
        interactive = true;
    end
end

if interactive,
    hGui = findobj('Tag','ObsviewGui');   % get handle to obsview gui, if it exists
else
    hGui = [];
end

% possibly remove data that failed QC

if ~isempty(hGui)&&strcmp(get(findobj(hGui,'Tag','qcxclear'),'Checked'),'on'),
    ods = odssubset(ods,ods.qcx==0);
    ods.ssqcx = 0;
end

nall = length(ods.kt);
if nall==0, return, end   % must have data to plot.

% decide whether to create a new figure window or plot in an existing window

if ~isempty(hGui), newfig = strcmp(get(findobj(hGui,'Tag','newfig'),'Checked'),'on'); 
else, newfig = 0; end

if ~newfig,         % see if we can find a previous obsplot figure window
    hFig = gcbf;
    if ~strcmp(get(hFig,'Tag'),'obsplot'), 
        h = findobj('Type','figure','Tag','obsplot'); 
        if any(h), hFig = h(end); else, newfig = 1; end
    end
end

if newfig,   
    hFig = figure;           % create a new figure window ...
    set(hFig,'Tag','obsplot','Color','w',...
        'DefaultTextFontUnits','normalized','DefaultTextFontSize',0.08)
    if interactive, 
        scr = get(0,'Screensize'); psize = min(0.9*scr(3:4)); 
        set(hFig,'Units','pixels','Position',psize*[0.1 0.1 1 8.5/11],...
            'MenuBar','none','DockControls','off')
        watchon;     % this might take a while
    end
else                         
    figure(hFig)             % ... or use an existing figure
    clf            
end

% set coloring indices

if ~isempty(hGui), 
    coloring = getappdata(hGui,'Coloring');
    preserve = strcmp(get(findobj(hGui,'Tag','keepcolors'),'Checked'),'on');
else   % defaults:
    if length(unique(ods.kx))==1
        coloring.attr = 'qcx';
    else  
        coloring.attr = 'kx';
    end
    coloring.cmap = 'jet';
    preserve = false;
end
if ~(isfield(ods,'cidx') & preserve), 
    ods = setcidx(ods,coloring); 
end

[KTPRS,KTSFC,KTRAD,SENSORS] = dconfig('KTPRS','KTSFC','KTRAD','SENSORS');

% count pressure-level/surface/radiance data

nprs = sum(ismember(ods.kt,KTPRS));
nsfc = sum(ismember(ods.kt,KTSFC));
nrad = sum(ismember(ods.kt,KTRAD));

% see which sensors are present in the data

sns = [];
for i = 1:length(SENSORS),
  if any(ismember(SENSORS(i).kxs,ods.kx)), sns = [sns i]; end
end

% create the different components of the figure

TtlPosition = [0.05 0.92 0.90 0.05];   % titles
MapPosition = [0.05 0.45 0.65 0.40];   % map with obs locations
LevPosition = [0.73 0.45 0.16 0.40];   % pressure distribution 
ChnPosition = [0.73 0.45 0.16 0.40];   % channel distribution
KtsPosition = [0.05 0.25 0.40 0.12];   % data type distribution
KxsPosition = [0.05 0.05 0.40 0.12];   % data source distribution
QcxPosition = [0.50 0.25 0.20 0.12];   % QC exclusion flag distribution
QchPosition = [0.73 0.25 0.20 0.12];   % QC history flag distribution
SnsPosition = [0.50 0.05 0.40 0.12];   % radiance data

hMapAxes = axes('Parent',hFig,'Tag','map','Position',MapPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.04);
figure(hFig), title('Observation Locations','Color','b','FontSize',0.06);
pmap('create',hMapAxes,ods)

if nprs>0,               % ods contains pressure-level data
    hLevAxes = axes('Parent',hFig,'Tag','lev','Position',LevPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.04);
    figure(hFig)
    title({'Vertical Distribution','of Pressure-Level Data'},...
        'Color','b','FontSize',0.06);
    plev('create',hLevAxes,ods)
elseif nrad>0 & length(sns)==1,   % radiance data from a single sensor
    hChnAxes = axes('Parent',hFig,'Tag','chn','Position',ChnPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.04);
    figure(hFig)
    title({'Channel Distribution',['for ' SENSORS(sns).id ' Data']},...
        'Color','b','FontSize',0.06);
    pchn('create',hChnAxes,ods)
end

hKtsAxes = axes('Parent',hFig,'Tag','kt','Position',KtsPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.12);
figure(hFig), title('Data Types','Color','b','FontSize',0.18)
pkts('create',hKtsAxes,ods)

hKxsAxes = axes('Parent',hFig,'Tag','kx','Position',KxsPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.12);
figure(hFig), title('Data Sources','Color','b','FontSize',0.18)
pkxs('create',hKxsAxes,ods)

if nrad>0,               % ods contains radiance data
    hSnsAxes = axes('Parent',hFig,'Tag','sn','Position',SnsPosition,'NextPlot','add',...
        'FontUnits','normalized','FontSize',0.12);
    figure(hFig), title('Sensor Distribution for Radiance Data',...
        'Color','b','FontSize',0.18)
    psns('create',hSnsAxes,ods)
end

hQcxAxes = axes('Parent',hFig,'Tag','qcx','Position',QcxPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.12);
figure(hFig), title('QC Exclusion','Color','b','FontSize',0.18)
pqcx('create',hQcxAxes,ods)

if ~ishandle(hQcxAxes), QchPosition = [0.50 0.25 0.40 0.12]; end
hQchAxes = axes('Parent',hFig,'Tag','qch','Position',QchPosition,'NextPlot','add',...
                'FontUnits','normalized','FontSize',0.12);
figure(hFig), title('QC History','Color','b','FontSize',0.18)
pqch('create',hQchAxes,ods)

hTtlAxes = axes('Parent',hFig,'Tag','title','Position',TtlPosition,'NextPlot','add',...
    'Visible','off',...
    'DefaultTextFontUnits','normalized','DefaultTextFontSize',0.4);
ptitle(hTtlAxes,ods,nall,nprs,nsfc,nrad)

if interactive
    
    % define some mouse actions

    set(hFig,'WindowButtonMotionFcn','obsplot pmap pointer')
    set([hMapAxes; allchild(hMapAxes)],'ButtonDownFcn','obsmap pobslist')

    % store the data with the plot

    set(hFig,'Userdata',ods)

    % create toolbar and menubar

    mktoolbar(hFig)
    mkmenubar(hFig)

    % create color legend

    plegend('create',hMapAxes,ods.cinfo,'NorthEast')

    watchoff;                            % done

end

%-------------------------------------------------------------------------
function panel(action,varargin)

switch action

case 'newfig'
    
    tag = varargin{1};

    hFig = gcbf;
    hAxes = findobj(hFig,'Type','axes','Tag',tag);

    hf = figure('Units','normalized','MenuBar','none',...
        'DockControls','off',...
        'Position',[0.1 0.2 0.8 0.8],'Visible','off');
    mktoolbar(hf)
    mkmenubar(hf)
    delete(findobj(hf,'Type','uimenu','Tag','Copy'))
    delete(findobj(hf,'Type','uimenu','Tag','Delete'))
    
    ha = copyobj(hAxes,hf);
    set(ha,'Units','normalized','Position',[0.1 0.15 0.8 0.7])
    ht = findobj(ha,'Type','text'); nt = numel(ht); % text objects
    fs = max(8,14-round(nt/10)); % adjust for crowded plots
    set(ht,'FontUnits','points','FontSize',fs);
    set(ha,'FontUnits','points','FontSize',fs); % tick labels
    set(findall(ha,'Type','text','HandleVisibility','off'),...
        'FontUnits','points','FontSize',20)  % title
    set(hf,'Userdata',get(hFig,'Userdata'))
    
    if strcmp(tag,'map'), % extra work for map panel

        hMenu = findobj(hFig,'Type','uimenu','Tag','showlegend');
        hNewMenu = findobj(hf,'Type','uimenu','Tag','showlegend');
        set(hNewMenu,'Checked',get(hMenu,'Checked'));
        ods = get(hf,'Userdata');
        plegend('create',ha,ods.cinfo,'NorthEastOutside')

        xlim = get(ha,'XLim'); ylim = get(ha,'YLim');
        asr = max(7/10,min(10/7,diff(xlim)*cos(mean(ylim)*pi/180)/diff(ylim)));
        w = 0.8; h = 0.8;
        if asr>1, h = h/asr; else, w = w*asr; end
        set(hf,'Position',[0.1 0.2 w h],...
            'WindowButtonMotionFcn','obsplot pmap pointer')
        set(hf,'Units','pixels')

    end

    set(hf,'Visible','on')
    
case 'delete'
            
    tag = varargin{1};
    
    hFig = gcbf;
    delete(findobj(hFig,'Type','axes','Tag',tag));
    set(findobj(hFig,'Type','uimenu','Tag',tag),'enable','off')

end

%-------------------------------------------------------------------------
function pmap(action,varargin)

switch action

case 'create'

    hAxes  = varargin{1};
    ods    = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)

    if isfield(ods,'sslon')&&any(ods.sslon), xlim = ods.sslon; else, xlim = [-180 180]; end
    if isfield(ods,'sslat')&&any(ods.sslat), ylim = ods.sslat; else, ylim = [ -90  90]; end
    
    % determine longitude limits xlim such that:
    %                -360<xlim(1)<xlim(2)<360 and xlim(2)-xlim(1)<=360
   
    xlim = -180 + mod(xlim + 180,360);     % xlim in [-180,180)
    
    if xlim(1)==xlim(2),                   % this means all longitudes: global map
        clon = 0;                          % by default, center global map at lon=0
        hGui = findobj('Tag','ObsviewGui');
        if any(hGui),                      % possible override by obsview gui option
            if strcmp(get(findobj(hGui,'Tag','dateline'),'Checked'),'on'), clon = 180; end
        end
        xlim = clon + [-180 180]; 
    elseif 0<=xlim(2)&&xlim(2)<xlim(1),
        xlim(1) = xlim(1) - 360;     
    elseif xlim(2)<0 &&xlim(2)<xlim(1),
        xlim(2) = xlim(2) + 360;         
    end                           
    
    set(hAxes,'XLim',xlim,'YLim',ylim,...
        'Box','on','XGrid','on','YGrid','on',...
        'Color',[0.8 0.9 1.0]);
    if diff(xlim)==360, set(hAxes,'XTick',linspace(xlim(1),xlim(2),7)); end
    if diff(ylim)==180, set(hAxes,'YTick',-60:30:60); end
    
    pcoast(hAxes)
    
    lon = ods.lon;
    if xlim(1)<-180, i = (lon>xlim(1)+360); lon(i) = lon(i) - 360; end
    if xlim(2)> 180, i = (lon<xlim(2)-360); lon(i) = lon(i) + 360; end

    iz = 0;
    ncolors = length(ods.cinfo);
    cvec = [3 2 1];
    if ncolors==3,
      for ax=1:ncolors,
        ix = cvec(ax);
        ic = ods.cidx==ix;
        if any(ic),
            iz = iz + 1;
            hLine = plot(lon(ic), ods.lat(ic), ...
                '.','Markersize',12,'Color',ods.cinfo(ix).rgb,'Tag','obsloc');
            set(hLine,'UserData',ix,'ZData',iz*ones(1,sum(ic)))
        end
      end
    else
      for ix = 1:length(ods.cinfo),
        ic = ods.cidx==ix;
        if any(ic),
            iz = iz + 1;
            hLine = plot(lon(ic), ods.lat(ic), ...
                '.','Markersize',12,'Color',ods.cinfo(ix).rgb,'Tag','obsloc');
            set(hLine,'UserData',ix,'ZData',iz*ones(1,sum(ic)))
%           disp([int2str(iz) ':  ' ods.cinfo(ix).txt])
        end
      end
    end
    
case 'domain'

    [h,hFig] = gcbo;                    % handles to cbo (h) and figure (hFig)
    hAxes = findobj(hFig,'Type','axes','Tag','map');  % handle to map axes

    region = lower(varargin{1});

    switch region
    case {'zoom'},               
        xlim = get(hAxes,'XLim');  ylim = get(hAxes,'YLim');
        if isequal(xlim,[-180 180])&&isequal(ylim,[-90 90]),
            warndlg({'You must first zoom in on the region you want:' ...
                     ' ' ...
                     'Press the zoom button on the figure window toolbar.' ...
                     'Then click, or click and drag, on the map to zoom.' ...
                     'Press the button again to get out of zoom mode.' ...
                     ' ' ...
                     'Then select this option again.'}, ...
                     'How to select data by zoom')
            return
        end
        xlim = -180 + mod(xlim + 180,360);   % xlim in [-180,180)
        if xlim(1)==xlim(2), xlim = [-180 180]; end
    case {'northernhemisphere(>20n)'}, xlim = [-180  180]; ylim = [  20  90];
    case {'southernhemisphere(>20s)'}, xlim = [-180  180]; ylim = [ -90 -20];
    case {'tropics(20s-20n)'},         xlim = [-180  180]; ylim = [ -20  20];
    case {'northamerica'},       xlim = [-172  -52]; ylim = [  16  72];
    case {'unitedstates'},       xlim = [-124  -72]; ylim = [  25  50];
    case {'northatlantic'},      xlim = [ -70    0]; ylim = [  10  70];
    case {'europe'},             xlim = [ -13   25]; ylim = [  35  72];
    case {'asia'},               xlim = [  25  170]; ylim = [   5  80];
    case {'arctica'},            xlim = [-180  180]; ylim = [  60  90];
    case {'pacific'},            xlim = [ 130 -120]; ylim = [ -60  40];
    case {'easternpacific'},     xlim = [ 130  180]; ylim = [   0  40];
    case {'westernpacific'},     xlim = [-180 -120]; ylim = [ -60  40];
    case {'southamerica'},       xlim = [ -85  -30]; ylim = [ -60  15];
    case {'southatlantic'},      xlim = [ -50   10]; ylim = [ -70   0];
    case {'africa'},             xlim = [ -20   52]; ylim = [ -40  38];
    case {'indianocean'},        xlim = [  40  100]; ylim = [ -40  10];
    case {'australia'},          xlim = [ 110  160]; ylim = [ -45 -10];
    case {'indonesia'},          xlim = [  90  130]; ylim = [ -10  10];
    case {'antarctica'},         xlim = [-180  180]; ylim = [ -90 -60];
    case {'amazon'},             xlim = [ -75  -45]; ylim = [ -12   6];
    case {'loweramazonbasin'},   xlim = [ -68  -57]; ylim = [ -22  -5];
    otherwise,                   xlim = [-180  180]; ylim = [ -90  90];

    end

    ods = get(hFig,'Userdata');            % get ods stored with the map
    ods.sslon = xlim;
    ods.sslat = ylim;
    if xlim(1)<xlim(2),
        odss = odssubset(ods, xlim(1)<=ods.lon&ods.lon<=xlim(2) &ylim(1)<=ods.lat&ods.lat<=ylim(2));
    else
        odss = odssubset(ods,(xlim(1)<=ods.lon|ods.lon<=xlim(2))&ylim(1)<=ods.lat&ods.lat<=ylim(2));
    end
    obsplot(odss)
        
case 'pointer'
    
    [h,hFig] = gcbo;                       % handles to cbo (h) and figure (hFig)

    % get current pointer position in normalized units:

    xy = get(hFig,'CurrentPoint');
    p = get(hFig,'Position');
    x = xy(1)/p(3);
    y = xy(2)/p(4);

    % see if the pointer is over the legend, if one is visible:

    hLeg = findobj(hFig,'Type','axes','Tag','legend');
    if ishandle(hLeg) & strcmp(get(hLeg,'Visible'),'on')
        p = get(hLeg,'Position');
        xl = (x-p(1))/p(3);
        yl = (y-p(2))/p(4);
        if (0<xl && xl<1 && 0<yl && yl<1),
           hLegTxt = findobj(hLeg,'Type','text');
           np = length(hLegTxt);
           ip = min(np,max(1,round(np*yl)));
           set(hLegTxt,'color','k')
           set(hLegTxt(ip),'color','r')
        end
    end

    % see if the pointer is over the map:

    hAxes = findobj(hFig,'Type','axes','Tag','map');     % handle to map axes
    hPloc = findobj(hAxes,'Type','text','Tag','pcoord'); % text object showing pointer position
    p = get(hAxes,'Position');
    xm = (x-p(1))/p(3);
    ym = (y-p(2))/p(4);
    if (0<xm && xm<1 && 0<ym && ym<1),

        % get lat-lon coordinates at pointer location

        xlim = get(hAxes,'XLim'); lon = xlim(1) + xm*(xlim(2)-xlim(1));
        ylim = get(hAxes,'YLim'); lat = ylim(1) + ym*(ylim(2)-ylim(1));
        str = [lonstr(lon) ', ' latstr(lat)];

        % modify or create text object showing pointer position

        if ishandle(hPloc),
            set(hPloc,'string',str)
        else
            text(xlim(1),ylim(2),str,'Parent',hAxes,'Tag','pcoord',...
                'HorizontalAlignment','left','VerticalAlignment','bottom',...
                'Color','r','Fontsize',0.04);
        end

        set(hFig,'Pointer','crosshair')

    else

        delete(hPloc)
        set(hFig,'Pointer','arrow')

    end

end
%-------------------------------------------------------------------------

function plegend(action,varargin)

switch action

case 'create'
    
    hAxes = varargin{1};
    cinfo = varargin{2};
    loc   = varargin{3};

    hLine = findobj(hAxes,'Type','line','Tag','obsloc'); % want legends for these objects
    nl = length(hLine);
    for il = 1:nl  
        ix = get(hLine(il),'UserData');    % color index
        legtxt{il} = cinfo(ix).txt;   % legend for this color
    end
    
    hLeg = legend(hAxes,hLine,legtxt,'Location',loc); % create the legend axes
    set(hLeg,'Color','w','UiContextMenu',[])
    
    hLegTxt = findobj(hLeg,'Type','text');  % context menu for each legend
    np = length(hLegTxt);
    for ip = 1:np
        ipstr = int2str(np-ip+1);
        visible = strcmp(get(hLine(np-ip+1),'Visible'),'on');
        hmenu = uicontextmenu;
        uimenu(hmenu,'Label','Show only','Tag','showonly',...
            'Callback',['obsplot plegend showonly ' ipstr])
        if visible, label = 'Hide'; else, label = 'Show'; end
        uimenu(hmenu,'Label',label,'Tag','showorhide',...
            'Callback',['obsplot plegend showorhide ' ipstr])
        uimenu(hmenu,'Separator','on','Label','Move to top',...
            'Callback',['obsplot plegend top ' ipstr])
        uimenu(hmenu,'Label','Move to bottom',...
            'Callback',['obsplot plegend bottom ' ipstr])
        uimenu(hmenu,'Separator','on','Label','Show all',...
            'Callback','obsplot plegend showorhide all')
        uimenu(hmenu,'Label','Hide all',...
            'Callback','obsplot plegend showorhide all')
        set(hLegTxt(ip),'UiContextMenu',hmenu,'UserData',hLine(ip))
    end
    
    hMenu = findobj(get(hLeg,'Parent'),'Type','uimenu','Tag','showlegend');  % handle to menu item
    on = strcmp(get(hMenu,'Checked'),'on');
    if ~on, legend(hAxes,'hide'); end
    
case 'showonly'
    
    ip = eval(varargin{1});
    
    [h,hFig] = gcbo;
    hLeg = findobj(hFig,'Type','axes','Tag','legend'); % handle to legend axes    
    hLegTxt = findobj(hLeg,'Type','text');  % legend text handles
    
    np = length(hLegTxt);
    for i = 1:np
        hp = get(hLegTxt(i),'UserData');  % line to be shown or hidden 
        hm = findobj(get(hLegTxt(np-i+1),'UiContextMenu'),'Tag','showorhide');
        if i==ip,
            set(hp,'Visible','on')
            set(hm,'Label','Hide')
        else
            set(hp,'Visible','off')
            set(hm,'Label','Show')
        end
    end  
    
case 'showorhide'
    
    [h,hFig] = gcbo;
    hLeg = findobj(hFig,'Type','axes','Tag','legend'); % handle to legend axes    
    hLegTxt = findobj(hLeg,'Type','text');  % legend text handles
    
    np = length(hLegTxt);
    if strcmp(varargin{1},'all')
        ip = 1:np;  % all items
    else
        ip = eval(varargin{1});  % single item
    end
    
    show = strncmp(get(h,'Label'),'Show',4);
    for i = ip
        hp = get(hLegTxt(i),'UserData');  % line to be shown or hidden 
        hm = findobj(get(hLegTxt(np-i+1),'UiContextMenu'),'Tag','showorhide');
        if show,
            set(hp,'Visible','on')
            set(hm,'Label','Hide')
        else
            set(hp,'Visible','off')
            set(hm,'Label','Show')
        end
    end

case 'showlegend'
    
    [h,hFig] = gcbo;
    hAxes = findobj(hFig,'Type','axes','Tag','map'); % handle to map axes
    hMenu = findobj(hFig,'Type','uimenu','Tag','showlegend');  % handle to menu item
    show = ~strcmp(get(hMenu,'Checked'),'on');
    if show
        legend(hAxes,'show')
        set(hMenu,'Checked','on')
    else
        legend(hAxes,'hide')
        set(hMenu,'Checked','off')
    end

case {'top','bottom'}
    
    ip = eval(varargin{1});  % sequence number of item to be moved
    
    [h,hFig] = gcbo;
    hAxes = findobj(hFig,'Type','axes','Tag','map');   % handle to map axes   
        
    hLeg = findobj(hFig,'Type','axes','Tag','legend'); % handle to legend axes    
    hLegTxt = findobj(hLeg,'Type','text');  % legend text handles 
  
    hp = get(hLegTxt(ip),'UserData');  % line to be moved
    zp = get(hp,'ZData'); zp = zp(1);
    np = length(hLegTxt); 
    
    if strcmp(action,'top'),
       for i = 1:np  
           h = get(hLegTxt(i),'UserData');
           z = get(h,'ZData');
           if z(1)>zp, set(h,'ZData',z-1); end
       end
       set(hp,'ZData',np*ones(size(get(hp,'ZData'))))
    else
       for i = 1:np  
           h = get(hLegTxt(i),'UserData');
           z = get(h,'ZData');
           if z(1)<zp, set(h,'ZData',z+1); end
       end
       set(hp,'ZData',ones(size(get(hp,'ZData'))))
    end   
    
%     ods = get(hFig,'UserData');
%     for i = 1:np
%         h = get(hLegTxt(i),'UserData');
%         z = get(h,'ZData');
%         ix = get(h,'UserData');
%         disp([int2str(z(1)) ':  ' ods.cinfo(ix).txt])
%     end

end

%-------------------------------------------------------------------------
function pcoast(hAxes)

xlim = get(hAxes,'XLim');

load coast
i = [0; find(isnan(long)); length(long)];
for j = 2:length(i)-1,
    is = i(j)+1:i(j+1)-1;
    patch(long(is),lat(is),[0.6 0.7 0.5])
end
for j = 1,     % this is a hack to fix Antartica
    is = i(j)+1:i(j+1)-1;
    patch([-180; long(is); 180],[-90; lat(is); -90],[0.6 0.7 0.5])
end

if xlim(1)<-180,  % this is a crude way to handle dateline crossing
    long = long - 360; 
    for j = 2:length(i)-1,
        is = i(j)+1:i(j+1)-1;
        patch(long(is),lat(is),[0.6 0.7 0.5])
    end
    for j = 1,
        is = i(j)+1:i(j+1)-1;
        patch([xlim(1); long(is); -180],[-90; lat(is); -90],[0.6 0.7 0.5])
    end
end

if xlim(2)>180,  % this is a crude way to handle dateline crossing
    long = long + 360; 
    for j = 2:length(i)-1,
        is = i(j)+1:i(j+1)-1;
        patch(long(is),lat(is),[0.6 0.7 0.5])
    end
    for j = 1,
        is = i(j)+1:i(j+1)-1;
        patch([180; long(is); xlim(2)],[-90; lat(is); -90],[0.6 0.7 0.5])
    end
end

%-------------------------------------------------------------------------
function plev(action,varargin)

switch action

case 'create'    

    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)
    
    KTPRS = dconfig('KTPRS');
    ipr = ismember(ods.kt,KTPRS);
    
    npmax = 35;
    p  = sort(unique(ods.lev(ipr)));
    np = length(p);

    if np<=1, delete(hAxes); return; end   % single level; forget it.

    logp = log(ods.lev(ipr));
    cidx = ods.cidx(ipr);

    % define bin centers z

    if np<npmax,       % discrete, irregularly spaced bins (in log(p))
        z = log(p);                         % bin centers
        pbins = [p' Inf];                   % bin edges [ )
    else,              % continuous, uniformly spaced bins (in log(p))
        [n,z] = hist(logp,npmax);           % bin centers
        maxlogp = max(logp); minlogp = min(logp);
        dz = (maxlogp-minlogp)/npmax;
        zz = minlogp + dz*(0:npmax);
        zz(end) = maxlogp;
        zbins = zz + max(eps,eps*abs(zz));
        pbins = exp(zbins);
        pbins(1) = 0; pbins(end) = Inf;     % bin edges [ )
    end
    
    nz = length(z);

    % count, by color
    
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nz);
    for ix = 1:nc,
       ncidx(ix,:) = hist(logp(cidx==ix),z);
       crgb{ix} = cinfo(ix).rgb;
    end

    % plot
    
    hp = barh(ncidx','stack');
    set(hp,{'FaceColor'},crgb')
    set(hAxes,'Box','on','YDir','reverse','YLim',[0 nz+1])

    % ticks

    ytick = 1:nz;
    if nz>20, ytick = 1:2:nz; end
    yticklabel = [];
    if np<npmax,
        for i = 1:length(ytick),
            yticklabel = strvcat(yticklabel,[num2str(exp(z(ytick(i))),4) 'hPa']);
        end
    else
        for i = 1:length(ytick),
            yticklabel = strvcat(yticklabel,[num2str(round(exp(z(ytick(i))))) 'hPa']);
        end
    end

    set(hAxes,'YTick',ytick,'YTickLabel',yticklabel,'YAxisLocation','right')
    set(hAxes,'XTick',[])
    
    % prepare callback

    set(hAxes,'Userdata',pbins)    % store bin edges with axes object
    set(hp,'ButtonDownFcn','obsplot plev barclick')

case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    pbins = get(hAxes,'Userdata'); % edges of pressure bins
    i = max(1,min(round(y),length(pbins)-1)); % bar index
       
    KTPRS = dconfig('KTPRS');
    s.kt  = KTPRS(ismember(KTPRS,ods.kt));
    s.lev = pbins([i i+1]);
    odss  = odssubset(ods,s);
    odss.sskt  = s.kt;
    minlev = min(odss.lev); maxlev = max(odss.lev);
    if diff([minlev maxlev]), 
        odss.sslev = s.lev;        % selected bin edges
    else
        odss.sslev = minlev;       % single level
    end  
    obsplot(odss)

end

%-------------------------------------------------------------------------
function pchn(action,varargin)

switch action

case 'create'
    
    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)
    
    ncmax = 35;
    ch = sort(unique(ods.lev));
    nc = length(ch);

    if nc<=1, delete(hAxes); return; end   % single channel; forget it.

    % define bin centers z

    if nc<ncmax,       % discrete, irregularly spaced bins
        z = ch;                         % bin centers
        chbins = [ch' Inf];             % bin edges [ )
    else,              % continuous, uniformly spaced bins
        [n,z] = hist(ch,ncmax);         % bin centers
        maxch = max(ch); minch = min(ch);
        dz = (maxch-minch)/ncmax;
        zz = minch + dz*(0:ncmax);
        zz(end) = maxch;
        chbins = zz;
        chbins(1) = 0; chbins(end) = Inf;     % bin edges [ )
    end
    
    nz = length(z);

    % count, by color
    
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nz);
    for ix = 1:nc,
       ncidx(ix,:) = hist(ods.lev(ods.cidx==ix),z);
       crgb{ix} = cinfo(ix).rgb;
    end

    % plot
    
    hp = barh(ncidx','stack');
    set(hp,{'FaceColor'},crgb')
    set(hAxes,'Box','on','YLim',[0 nz+1])

    % ticks

    ytick = 1:nz;
    if nz>20, ytick = 1:2:nz; end
    yticklabel = [];
    if nc<ncmax,
        for i = 1:length(ytick),
            yticklabel = strvcat(yticklabel,num2str(z(ytick(i))));
        end
    else
        for i = 1:length(ytick),
            yticklabel = strvcat(yticklabel,num2str(round(z(ytick(i)))));
        end
    end

    set(hAxes,'YTick',ytick,'YTickLabel',yticklabel,'YAxisLocation','right')
    set(hAxes,'XTick',[])
    
    % prepare callback

    set(hAxes,'Userdata',chbins)    % store bin edges with axes object
    set(hp,'ButtonDownFcn','obsplot pchn barclick')

case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    pbins = get(hAxes,'Userdata'); % edges of pressure bins
    i = max(1,min(round(y),length(pbins)-1)); % bar index

    odss = odssubset(ods,pbins(i)<=ods.lev & ods.lev<pbins(i+1));
    minlev = min(odss.lev); maxlev = max(odss.lev);
    if diff([minlev maxlev]), 
        odss.sslev = pbins([i i+1]); 
    else, 
        odss.sslev = minlev;       % single channel
    end
    obsplot(odss)

end

%-------------------------------------------------------------------------
function pkts(action,varargin)

switch action

case 'create'

    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)

    kts = unique(ods.kt);
    nkt = length(kts);
    
    % count, by color
    
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nkt);
    if nkt==0,
        delete(hAxes), return
    elseif nkt==1,
        for ix = 1:nc,
            ncidx(ix,:) = sum(ods.cidx==ix);
            crgb{ix} = cinfo(ix).rgb;
        end
    else
        for ix = 1:nc,
            ncidx(ix,:) = hist(single(ods.kt(ods.cidx==ix)),single(kts));
            crgb{ix} = cinfo(ix).rgb;
        end
    end
    n = sum(ncidx,1);

    % plot
    
    hp = barh([ncidx';zeros(1,nc)],'stack');
    set(hp,{'FaceColor'},crgb')  
    hcb = hp;

    set(hAxes,'YLim',[0 nkt+1],'YTick',[])

    xtick = [0.50 1.00 1.50]*max(n);
    factor = max(1,10^floor(log10(max(n))));
    for i = 1:length(xtick),
        xt = xtick(i)/factor;
        xticklabel{i} = num2str(xt,2);
    end
    if factor>1,
        xticklabel{end} = ['(x ' int2str(factor) ')'];
    end
    set(hAxes,'XLim',[0 1.5*max(n)],'XTick',xtick,'XTickLabel',xticklabel)

    KTS = dconfig('KTS');
    fs = min([0.17 max([0.07 0.7/nkt])]);
    for i = 1:nkt,
        j = find(kts(i)==[KTS.value]);
        if any(j),
            hcb = [hcb text(n(i),i,['  ' KTS(j).id],'FontSize',fs)];
        end
    end
    set(hAxes,'YTick',1:nkt,'YTickLabel',kts)
    
    % prepare callback

    if nkt>1,
        set(hAxes,'Userdata',kts)  % store kts with axes object
        set(hcb,'ButtonDownFcn','obsplot pkts barclick'); 
    end   
    
case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    kts = get(hAxes,'Userdata');   % kts associated with bars
    i = max(1,min(round(y),length(kts))); % bar index

    ods.sskt = kts(i);
    obsplot(odssubset(ods,ods.kt==kts(i)))

end

%-------------------------------------------------------------------------
function pkxs(action,varargin)

switch action

case 'create'

    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)

    kxs = unique(ods.kx);
    nkx = length(kxs);

    % count, by color
    
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nkx);
    if nkx==0,
        delete(hAxes), return
    elseif nkx==1,
        for ix = 1:nc,
            ncidx(ix,:) = sum(ods.cidx==ix);
            crgb{ix} = cinfo(ix).rgb;
        end
    else
        for ix = 1:nc,
            ncidx(ix,:) = hist(single(ods.kx(ods.cidx==ix)),single(kxs));
            crgb{ix} = cinfo(ix).rgb;
        end
    end
    n = sum(ncidx,1);

    % plot
    
    hp = barh([ncidx';zeros(1,nc)],'stack');
    set(hp,{'FaceColor'},crgb')  
    hcb = hp;
    set(hAxes,'YLim',[0 nkx+1],'YTick',[])

    xtick = [0.50 1.00 1.50]*max(n);
    factor = max(1,10^floor(log10(max(n))));
    for i = 1:length(xtick),
        xt = xtick(i)/factor;
        xticklabel{i} = num2str(xt,2);
    end
    if factor>1,
        xticklabel{end} = ['(x ' int2str(factor) ')'];
    end
    set(hAxes,'XLim',[0 1.5*max(n)],'XTick',xtick,'XTickLabel',xticklabel)

    KXS = dconfig('KXS');
    fs = min([0.17 max([0.07 0.7/nkx])]);
    for i = 1:nkx,
        j = find(kxs(i)==[KXS.value]);
        if any(j),
            hcb = [hcb text(n(i),i,['  ' KXS(j).id],'FontSize',fs)];
        end
    end
    set(hAxes,'YTick',1:nkx,'YTickLabel',kxs)
    
    % prepare callback

    if nkx>1,
        set(hAxes,'Userdata',kxs)  % store kxs with axes object
        set(hcb,'ButtonDownFcn','obsplot pkxs barclick'); 
    end

case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    kxs = get(hAxes,'Userdata');   % kxs associated with bars
    i = max(1,min(round(y),length(kxs))); % bar index

    ods.sskx = kxs(i);
    obsplot(odssubset(ods,ods.kx==kxs(i)))

end

%-------------------------------------------------------------------------
function psns(action,varargin)

switch action

case 'create'

    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)
    
    SENSORS = dconfig('SENSORS');
    
    sn = NaN*ones(size(ods.kx));
    for i = 1:length(SENSORS),
       sn(ismember(ods.kx,SENSORS(i).kxs)) = i;
    end
    isn = isfinite(sn);     % mask for data with identifiable sensors
    sns = unique(sn(isn));  % list of sensors present in the data
    nsn = length(sns);      % number of sensors present in the data

    sn   = sn(isn);
    cidx = ods.cidx(isn);
        
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nsn);
    if nsn==0,
        delete(hAxes), return
    elseif nsn==1,
        for ix = 1:nc,
            ncidx(ix,:) = sum(cidx==ix);
            crgb{ix} = cinfo(ix).rgb;
        end
    else
        for ix = 1:nc,
            ncidx(ix,:) = hist(single(sn(cidx==ix)),single(sns));
            crgb{ix} = cinfo(ix).rgb;
        end
    end
    n = sum(ncidx,1);
    
    % plot
    
    hp = barh([ncidx';zeros(1,nc)],'stack');
    set(hp,{'FaceColor'},crgb')  
    hcb = hp;
    set(hAxes,'YLim',[0 nsn+1],'YTick',[])

    xtick = [0.50 1.00 1.50]*max(n);
    factor = max(1,10^floor(log10(max(n))));
    for i = 1:length(xtick),
        xt = xtick(i)/factor;
        xticklabel{i} = num2str(xt,2);
    end
    if factor>1,
        xticklabel{end} = ['(x ' int2str(factor) ')'];
    end
    set(hAxes,'XLim',[0 1.5*max(n)],'XTick',xtick,'XTickLabel',xticklabel)

    fs = min([0.17 max([0.07 0.7/nsn])]);
    for i = 1:nsn,
        hcb = [hcb text(n(i),i,['  ' SENSORS(sns(i)).id],'FontSize',fs)];
    end
    set(hAxes,'YTick',1:nsn,'YTickLabel',[])
    
    % prepare callback

    if nsn>1,
        set(hAxes,'Userdata',sns)  % store sns with axes object
        set(hcb,'ButtonDownFcn','obsplot psns barclick'); 
    end

case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    sns = get(hAxes,'Userdata');   % sensors associated with bars
    i = max(1,min(round(y),length(sns))); % bar index
    
    SENSORS = dconfig('SENSORS');
    s.kx = SENSORS(sns(i)).kxs;
    ods.sskx = s.kx;
    obsplot(odssubset(ods,s))

end

%-------------------------------------------------------------------------
function pqcx(action,varargin)

switch action

case 'create'

    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)

    iok  = ods.qcx>0; % consider nonzero qcx only
    qcx  = ods.qcx(iok); 
    cidx = ods.cidx(iok);
    qcs  = unique(qcx); 
    nqc  = length(qcs);
    
    % count, by color
    
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nqc);
    if nqc==0 % don't bother with this panel if there are no non-zero qc flags
        delete(hAxes), return
    elseif nqc==1,
        for ix = 1:nc,
            ncidx(ix,:) = sum(cidx==ix);
            crgb{ix} = cinfo(ix).rgb;
        end
    else
        for ix = 1:nc,
            ncidx(ix,:) = hist(single(qcx(cidx==ix)),single(qcs));
            crgb{ix} = cinfo(ix).rgb;
        end
    end
    n = sum(ncidx,1);

    % plot
    
    hp = barh([ncidx';zeros(1,nc)],'stack');
    set(hp,{'FaceColor'},crgb')  
    hcb = hp;
    set(hAxes,'YLim',[0 nqc+1],'YTick',[])

    xtick = [0.50 1.00 1.50]*max(n);
    factor = max(1,10^floor(log10(max(n))));
    for i = 1:length(xtick),
        xt = xtick(i)/factor;
        xticklabel{i} = num2str(xt,2);
    end
    if factor>1,
        xticklabel{end} = ['(x ' int2str(factor) ')'];
    end
    set(hAxes,'XLim',[0 1.5*max(n)],'XTick',xtick,'XTickLabel',xticklabel)

    QCXS = dconfig('QCXS');
    fs = min([0.17 max([0.07 0.7/nqc])]);
    for i = 1:length(QCXS),
        j = find(qcs==QCXS(i).value);
        if any(j),
            hcb = [hcb text(n(j),j,['  ' QCXS(i).id],'FontSize',fs)];
        end
    end
    set(hAxes,'YTick',1:nqc,'YTickLabel',qcs)
    
    % prepare callback

    set(hAxes,'Userdata',qcs)      % store qcs with axes object
    set(hcb,'ButtonDownFcn','obsplot pqcx barclick'); 

case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    qcs = get(hAxes,'Userdata');   % qcs associated with bars
    i = max(1,min(round(y),length(qcs))); % bar index

    ods.ssqcx = qcs(i);
    obsplot(odssubset(ods,ods.qcx==qcs(i)))

end

%-------------------------------------------------------------------------
function pqch(action,varargin)

switch action

case 'create'

    hAxes = varargin{1};
    ods   = varargin{2};

    hFig = get(hAxes,'Parent');
    figure(hFig)

    qch  = ods.qch; 
    cidx = ods.cidx;
    qcs  = unique(qch); 
    nqc  = length(qcs);
    
    % count, by color
    
    cinfo = ods.cinfo;
    nc = length(cinfo);
    ncidx = zeros(nc,nqc);
    if nqc==1 % don't bother with this panel if there is only 1 flag
        delete(hAxes), return
    else
        for ix = 1:nc,
            ncidx(ix,:) = hist(single(qch(cidx==ix)),single(qcs));
            crgb{ix} = cinfo(ix).rgb;
        end
    end
    n = sum(ncidx,1);

    % plot
    
    hp = barh([ncidx';zeros(1,nc)],'stack');
    set(hp,{'FaceColor'},crgb')  
    hcb = hp;
    set(hAxes,'YLim',[0 nqc+1],'YTick',[])

    xtick = [0.50 1.00 1.50]*max(n);
    factor = max(1,10^floor(log10(max(n))));
    for i = 1:length(xtick),
        xt = xtick(i)/factor;
        xticklabel{i} = num2str(xt,2);
    end
    if factor>1,
        xticklabel{end} = ['(x ' int2str(factor) ')'];
    end
    set(hAxes,'XLim',[0 1.5*max(n)],'XTick',xtick,'XTickLabel',xticklabel)

    QCHS = dconfig('QCHS');
    fs = min([0.17 max([0.07 0.7/nqc])]);
    for i = 1:length(QCHS),
        j = find(qcs==QCHS(i).value);
        if any(j),
            hcb = [hcb text(n(j),j,['  ' QCHS(i).id],'FontSize',fs)];
        end
    end
    set(hAxes,'YTick',1:nqc,'YTickLabel',qcs,'YDir','reverse')
    
    % prepare callback

    set(hAxes,'Userdata',qcs)      % store qcs with axes object
    set(hcb,'ButtonDownFcn','obsplot pqch barclick'); 

case 'barclick'

    [hp,hFig] = gcbo;              % handles to patch (hp) and figure (hFig)
    if strcmp(get(hFig,'SelectionType'),'alt'), return, end  % if right-click, do nothing
    hAxes = get(hp,'Parent');      % handle to axes
    v = get(hAxes,'CurrentPoint');
    y = v(1,2);                    % y-coordinate of mouse click

    ods = get(hFig,'Userdata');    % original data
    qcs = get(hAxes,'Userdata');   % qcs associated with bars
    i = max(1,min(round(y),length(qcs))); % bar index

    ods.ssqch = qcs(i);
    obsplot(odssubset(ods,ods.qch==qcs(i)))

end

%-------------------------------------------------------------------------
function ptitle(hAxes,ods,nall,nprs,nsfc,nrad)

ssid  = []; if isfield(ods,'ssid' ), ssid  = ods.ssid ; end
sslat = []; if isfield(ods,'sslat'), sslat = ods.sslat; end
sslon = []; if isfield(ods,'sslon'), sslon = ods.sslon; end
sslev = []; if isfield(ods,'sslev'), sslev = ods.sslev; end
sskt  = []; if isfield(ods,'sskt' ), sskt  = ods.sskt ; end
sskx  = []; if isfield(ods,'sskx' ), sskx  = ods.sskx ; end
ssqcx = []; if isfield(ods,'ssqcx'), ssqcx = ods.ssqcx; end
ssqch = []; if isfield(ods,'ssqch'), ssqch = ods.ssqch; end

line1 = [ssid ': ' int2str(nall) ' observations'];
if ~isempty(sslat),
    str = sprintf('lat=[%6.1f,%6.1f)',sslat);
    str = strrep(str,' ',''); str = strrep(str,'.0','');
else
    str = 'all lat';
end
line2 = [str ';'];
if ~isempty(sslon),
    str = sprintf('lon=[%6.1f,%6.1f)',sslon);
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
figure(hFig), axes(hAxes)
set(gca,'XLim',[0 1],'YLim',[0 1])
text(0,1,line1,'Tag','title1',...
    'Interpreter','none','Color','blue','HorizontalAlignment','left',...
    'FontSize',0.5)
text(0,0.33,line2,'Tag','title2',...
    'Interpreter','none','HorizontalAlignment','left')
if ~isempty(line3),
    text(0,0,line3,'Tag','title3',...
    'Interpreter','none','HorizontalAlignment','left');
end

% make a legend

text(1,1.00,sprintf('%s%10d','Radiance Data: ',nrad),...
        'Color','b','HorizontalAlignment','right','FontSize',0.3)
text(1,0.66,sprintf('%s%10d','Pressure-Level Data: ',nprs),...
        'Color','b','HorizontalAlignment','right','FontSize',0.3)
text(1,0.33,sprintf('%s%10d','Surface Data: ',nsfc),...
        'Color','b','HorizontalAlignment','right','FontSize',0.3)
noth = nall-nrad-nprs-nsfc;
text(1,0.00,sprintf('%s%10d','Other Data: ',noth),...
        'Color','b','HorizontalAlignment','right','FontSize',0.3)

%--------------------------------------------------------------------------
function select(attr)

[h,hFig] = gcbo;      % handles to callback object (h) and figure (hFig)
ods = get(hFig,'Userdata');   % get ods stored with the plot
title = 'Select by attribute';

[OBSATTRIBUTES,KTS,KXS,QCHS,QCXS] = ...
   dconfig('OBSATTRIBUTES','KTS','KXS','QCHS','QCXS');
attrs = OBSATTRIBUTES.names;
descr = OBSATTRIBUTES.descr;
discr = OBSATTRIBUTES.discr;

ia = find(strcmp(attr,attrs));
if any(ia)

   if discr(ia) % discrete attribute values

      values = unique(ods.(attr)); n = length(values);
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

      s.(attr) = values(i)'; % subset list of values

   else        % continuous attribute values

      minx = min(ods.(attr));
      maxx = max(ods.(attr));
      prompt  = {['Minimum ' descr{ia} ':'],['Maximum ' descr{ia} ':']};
      def     = {num2str(minx),num2str(maxx)};
      answer  = inputdlg(prompt,title,1,def);
      if isempty(answer), return; end
      minx    = eval(answer{1});
      maxx    = eval(answer{2});

      if minx==maxx, s.(attr) = minx;     % subset single value
      else, s.(attr) = sort([minx maxx]); % subset range of values
      end

   end

   ods.(['ss' attr]) = s.(attr);  % store subset criteria with ods
   obsplot(odssubset(ods,s))      % create the new plot

end

%--------------------------------------------------------------------------
function mktoolbar(hFig)

% Toolbar

set(0,'ShowHiddenHandles','on')
set(hFig,'Toolbar','figure') % need this for a moment
keep = {'Standard.SaveFigure';'Standard.PrintFigure';...
    'Exploration.ZoomIn';'Exploration.ZoomOut';'Exploration.Rotate'};
ht = uitoolbar;              % create custom toolbar:
for i = 1:length(keep)
    h = findobj(hFig,'Tag',keep{i});
    copyobj(h,ht);
end
set(0,'ShowHiddenHandles','off')
set(hFig,'Toolbar','none')   % remove the standard toolbar

%--------------------------------------------------------------------------
function mkmenubar(hFig)

% File menu 

h = uimenu(hFig,'Label','File'); 
uimenu(h,'Label','Print...',...
    'Callback','filemenufcn FilePrintSetup');
uimenu(h,'Label','Save...',...
    'Callback','filemenufcn FileSaveAs');
uimenu(h,'Label','Page Setup...',...
    'Callback','filemenufcn FilePageSetup');
uimenu(h,'Label','Close','Separator','on',...
    'Callback','filemenufcn FileClose');

% Subset menu

h = uimenu(hFig,'Label','Subset'); 
hs = uimenu(h,'Label','Based on attribute values...');
ATTRS = dconfig('OBSATTRIBUTES');
for i = 1:length(ATTRS.names),
    uimenu(hs,'Label',[ATTRS.names{i} ' (' ATTRS.descr{i} ')'],...
        'Callback',['obsplot select ' ATTRS.names{i}]);
end

hs = uimenu(h,'Label','Based on geographic region...');
for region = {'Northern Hemisphere (>20N)','Southern Hemisphere (>20S)','Tropics (20S-20N)',...
            'North America','United States','North Atlantic',...
            'Europe','Asia','Arctica',...
            'Pacific','Eastern Pacific','Western Pacific',...
            'South America','Amazon','Lower Amazon Basin',...
            'South Atlantic','Africa','Indian Ocean',...
            'Indonesia','Australia','Antarctica'},
    uimenu(hs,'Tag',region{:},'Label',region{:},...
        'Callback',['obsplot pmap domain ' strrep(region{:},' ','')]);
end
uimenu(h,'Label','Based on current zoom state',...
    'Callback','obsplot pmap domain zoom');
if isempty(findobj('Tag','ObsviewGui')), enable = 'off'; else, enable = 'on'; end
uimenu(h,'Label','Save subset definition...','Separator','on',...
    'Callback','obsview addgroup','Enable',enable);
uimenu(h,'Label','Send subset to workspace',...
    'Callback','ods=get(gcbf,''Userdata'')');

% Plot menu

h = uimenu(hFig,'Label','Plot'); 
hC = uimenu(h,'Tag','Copy','Label','Copy panel to a new figure...');
hD = uimenu(h,'Tag','Delete','Label','Delete panel from this figure...');
m(1).titl = 'Observation Locations'; m(1).tag = 'map'; m(1).pcall = 'pmap';
m(2).titl = 'Vertical Distribution'; m(2).tag = 'lev'; m(2).pcall = 'plev';
m(3).titl = 'Channel Distribution';  m(3).tag = 'chn'; m(3).pcall = 'pchn';
m(4).titl = 'Data Types';            m(4).tag = 'kt';  m(4).pcall = 'pkts';
m(5).titl = 'Data Sources';          m(5).tag = 'kx';  m(5).pcall = 'pkxs';
m(6).titl = 'Sensor Distribution';   m(6).tag = 'sn';  m(6).pcall = 'psns';
m(7).titl = 'QC Exclusion';          m(7).tag = 'qcx'; m(7).pcall = 'pqcx';
m(8).titl = 'QC History';            m(8).tag = 'qch'; m(8).pcall = 'pqch';
for i = 1:length(m),
    if isempty(findobj('Type','axes','Tag',m(i).tag)), enable = 'off'; else, enable = 'on'; end
    uimenu(hC,'Tag',m(i).tag,'Label',m(i).titl,'Enable',enable,...
                 'Callback',['obsplot panel newfig ' m(i).tag]);
    uimenu(hD,'Tag',m(i).tag,'Label',m(i).titl,'Enable',enable,...
                 'Callback',['obsplot panel delete ' m(i).tag]);
end
uimenu(h,'Label','Create interactive map','Separator','on',...
    'Callback','obsmap(get(gcbf,''Userdata''))');
uimenu(h,'Label','Create data list',...
    'Callback','obslist(get(gcbf,''Userdata''))')
uimenu(h,'Label','Create histogram...',...
    'Callback','obshist(get(gcbf,''Userdata''))')


% Options menu

hGui = findobj('Tag','ObsviewGui');   % get handle to obsview gui, if it exists

h = uimenu(hFig,'Label','Options'); 
if ~isempty(hGui), checked = get(findobj(hGui,'Tag','clegend'),'Checked'); else, checked = 'off'; end
uimenu(h,'Label','Show color legend',...
    'Tag','showlegend','Checked',checked,'Callback','obsplot plegend showlegend')

%------------------------------------------------------------------------
function toggle(option)

h = findobj(gcbf,'Type','uimenu','Tag',option);
if strcmp(get(h,'Checked'),'on'), set(h,'Checked','off');
else, set(h,'Checked','on'); end


