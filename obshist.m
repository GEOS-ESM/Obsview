function obshist(arg,varargin)

% OBSHIST Create a histogram of observation attributes.
%
% OBSHIST(ODS) brings up a list of attributes for the user to choose from,
%    and then produces a histogram for that attribute
% OBSHIST(ODS,ATTR) produces a histogram for the attribute ATTR

% 27Aug2002 Dick Dee (dee@dao.gsfc.nasa.gov)
% 23Apr2004 Dick Dee - ECMWF version

% check whether this is a callback or not

if ischar(arg), % must be a callback
   feval(arg,varargin{:})      % simply pass to the callback routine..
   return                      % .. and get out.
else            % otherwise this will set up a new plot
   ods = arg;  % first arg is always the ods structure with the data
end

% check if we have data to plot

n = length(ods.kt);
if n==0, return, end

OBSATTRIBUTES = dconfig('OBSATTRIBUTES');
attrs = OBSATTRIBUTES.names;
descr = OBSATTRIBUTES.descr;

% get the attribute name

attr = 'xxx';
if nargin>1, attr = varargin{1}; end

if ~strncmp(attr,attrs,length(attr)),     % ask the user
   for i = 1:length(attrs)
      choices{i} = [descr{i} ' (' attrs{i} ')'];
   end
   name   = 'Create a histogram';
   prompt = 'Select a data attribute:';
   i = listdlg('Name',name,'PromptString',prompt,...
      'ListSize',[200 200],'ListString',choices,'SelectionMode','single');
   if isempty(i), return; end
   attr = attrs{i};
end
c = double(ods.(attr));          % the data

% see if there is coloring info for these data

if isfield(ods,'cidx')&&isfield(ods,'cinfo'),
   cidx = ods.cidx;
   cinfo = ods.cinfo;
else
   cidx = ones(size(c),'int8');
   cinfo.rgb = [0 1 0];
   cinfo.txt = '';
end

junk = isnan(c);
if any(junk),
   c(junk) = [];                % remove junk
   cidx(junk) = [];
   if isempty(c), warndlg('No data for this attribute.'), return, end
   n = length(c);
end

x = unique(c);
cs = unique(cidx);
nc = length(cs); % number of colors

switch attr

   case {'kx','kt','qcx','qch'},
      continuous = 0;                                % discrete range
   otherwise
      continuous = 1;                                % continuous range
      if length(x)<20, continuous = 0; end           % (on second thought)

end

if continuous,

   xmin = min(c); xmax = max(c); nx = 50;
   x = linspace(xmin,xmax,nx);

   % count, by color

   ncidx = zeros(nc,nx);
   for ix = 1:nc,
      i = cs(ix);  % color index
      ncidx(ix,:) = hist(c(cidx==i),x);
      crgb{ix} = cinfo(i).rgb;
      legs{ix} = cinfo(i).txt;
   end
   ntotal = sum(ncidx,1);

   histFig = figure('Tag','obshist');
   figure(histFig)

   h = bar(x,ncidx','stack');
   set(h,'EdgeColor','k',{'FaceColor'},crgb')
   hold on

   % Create legend for this plot

   legend(h,legs);

   ddx = 0.05*(xmax-xmin);
   axis([xmin-ddx xmax+ddx 0 1.05*max(ntotal)])

   % Create (invisible) Gaussian

   fx = (n*(x(2)-x(1))/(std(c)*sqrt(2*pi)))*exp(-0.5*((x-mean(c))/std(c)).^2);
   plot(x,fx,'c-','LineWidth',2,'Visible','off','Tag','gaussian')
   hold off

   % Create a title for the plot

   title([attr ': mean=' num2str(mean(c)) ', stdv = ' num2str(std(c))],...
      'Color','b','Fontsize',16)

else

   if length(x)==1,
      x = x + (-1:1)';
   elseif length(x)>20,
      [n,x] = hist(c);
   end
   nx = length(x);

   % count, by color

   ncidx = zeros(nc,nx);
   for ix = 1:nc,
      i = cs(ix);  % color index
      ncidx(ix,:) = hist(c(cidx==i),x);
      crgb{ix} = cinfo(i).rgb;
      legs{ix} = cinfo(i).txt;
   end
   ntotal = sum(ncidx,1);

   histFig = figure('Tag','obshist');
   figure(histFig)

   h = bar(x,ncidx','stack');
   set(h,'EdgeColor','k',{'FaceColor'},crgb')
   hold on

   legend(h,legs);

   if length(x)>10,
      xmax = max(c);
      xmin = min(c);
      dx = 0.05*(xmax-xmin);
      set(gca,'XLim',[xmin-dx xmax+dx])
   end
   set(gca,'YLim',[0 1.05*max(ntotal)])
   title(attr,...
      'Color','b','Fontsize',16)

end

% Create context menus for the plot

hMenu = uicontextmenu;
set(gca,'UIContextMenu',hMenu);

uimenu(hMenu,'Label','Add or delete legend box',...
   'Callback','legend toggle');
if (continuous)
   uimenu(hMenu,'Label','Add a Gaussian curve',...
      'Callback','delete(gcbo); set(findobj(gca,''Tag'',''gaussian''),''Visible'',''on'')');
end
