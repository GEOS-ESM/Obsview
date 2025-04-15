function ptsbars(hfig,ts,ir,j,expid)

% PTSBARS Plot time series of statistics for OBS2HTML.

%     ts.id
%     ts.isprs
%     ts.israd
%     ts.reg  { :          }  
%     ts.lev  (     :      )
%     ts.time (         :  )
%     ts.nobs (reg,lev,time)  
%     ts.nana (reg,lev,time)
%     ts.npsv (reg,lev,time)
%     ts.somf (reg,lev,time)
%     ts.somf2(reg,lev,time)
%     ts.soma (reg,lev,time)
%     ts.soma2(reg,lev,time)
%     ts.Jobkg(reg,lev,time)
%     ts.Joana(reg,lev,time)

if nargin<4, expid = ''; end

% sort by time:

[time,is] = sort(ts.time);

nobs   = squeeze(ts.nobs (ir,j,is));
nana   = squeeze(ts.nana (ir,j,is));
npsv   = squeeze(ts.npsv (ir,j,is));
somf   = squeeze(ts.somf (ir,j,is));
somf2  = squeeze(ts.somf2(ir,j,is));
soma   = squeeze(ts.soma (ir,j,is));
soma2  = squeeze(ts.soma2(ir,j,is));
tJobkg = squeeze(ts.Jobkg(ir,j,is));
tJoana = squeeze(ts.Joana(ir,j,is));

nexcl = nobs - (nana + npsv);
nt = length(time);
zero = zeros(nt,1);

% create bar plots:

set(hfig,'Units','points','DefaultAxesFontSize',12)
figure(hfig), clf

% Plot title
% ----------
subplot('Position',[0.1 0.9 0.8 0.1])

   t1 = time(1); t2 = time(end);
   tstr = [datestr(t1,'dd') datestr(t1,'mmm') datestr(t1,'yyyy') ' ' datestr(t1,'HH') 'Z'];
   if t2>t1,
       tstr = [tstr ' - ' datestr(t2,'dd') datestr(t2,'mmm') datestr(t2,'yyyy') ' ' datestr(t2,'HH') 'Z'];
   end
   ptitle = tstr;
   if any(expid)
     ptitle = [expid '  ' ptitle];
   end
   text(0.5,0.5,ptitle,'FontSize',12,'FontWeight','b','Color','k',...
                       'HorizontalAlignment','center',...
                       'VerticalAlignment','bottom',...
                       'Interpreter','none');
   ptitle = ts.id;                
   if any(ts.lev)
      if     ts.isprs, levstr = [num2str(ts.lev(j)) 'hPa'];
      elseif ts.israd, levstr = ['channel ' int2str(ts.lev(j))];
      else             levstr = ['lev = ' int2str(ts.lev(j))]; end
      ptitle = [ptitle ':  ' levstr];
   end
   ptitle = [ptitle '  (' ts.reg{ir} ')'];
   text(0.5,0.1,ptitle,'FontSize',12,'FontWeight','b','Color','b',...
                       'HorizontalAlignment','center',...
                       'VerticalAlignment','bottom',...
                       'Interpreter','none');
   set(gca,'Visible','off')

% Data count
% ----------
subplot('Position',[0.1 0.65 0.8 0.22])

   if nt>1
       h = bar(time, [nana npsv nexcl],'stack');
   else
       h = bar([time time time; 0 0 0], [nana npsv nexcl; 0 0 0],'stack');
   end
   set(h,{'FaceColor'},{'g';[1 0.5 0];'r'}, 'EdgeColor','none')
   ylim = [0 1.1*max(nobs)];
   ptaxes(time,ylim,1,'Data counts:','k','Used (p)','g','Passive',[1 0.5 0],'Not used','r')

% Omf, oma statistics
% -------------------
subplot('Position',[0.1 0.35 0.8 0.22])

   n = nana + npsv;
   i = (n>0);
   trmsf = zero; trmsf(i) = sqrt(somf2(i)./n(i));
   trmsa = zero; trmsa(i) = sqrt(soma2(i)./n(i));
   tmeanf = zero; tmeanf(i) = somf(i)./n(i);
   tmeana = zero; tmeana(i) = soma(i)./n(i);
   nt = length(time);
   if nt>1, dt = 0.1*(time(end)-time(1))/(nt-1);
   else, dt = 0.2; end
   hrf = bar(time-dt,trmsf,0.6); hold on
   hra = bar(time+dt,trmsa,0.6); hold on
   hmf = bar(time-dt,tmeanf,0.6); hold on
   hma = bar(time+dt,tmeana,0.4); hold off
   set(hrf,'FaceColor','b','EdgeColor','none')
   set(hra,'FaceColor','r','EdgeColor','none')
   set(hmf,'FaceColor','c','EdgeColor','none')
   set(hma,'FaceColor',[1 0.8 0],'EdgeColor','none')
   line([time(1)-1 time(end)+1],[0 0],'Color','k')
   ylim = 1.1*[min([0 tmeanf' tmeana']) max([trmsf' trmsa'])];
   ptaxes(time,ylim,0,'Data residuals:','k','rms(O-B)','b','rms(O-A)','r','mean(O-B)','c','mean(O-A)',[1 0.8 0])
     
% Jo(bkg)/n, Jo(ana)/n
% --------------------
subplot('Position',[0.1 0.05 0.8 0.22])

   n = nana;
   i = (n>0);
   tJobkg(i) = tJobkg(i)./n(i);
   tJoana(i) = tJoana(i)./n(i);
   nt = length(time);
   if nt>1, dt = 0.1*(time(end)-time(1))/(nt-1);
   else, dt = 0.2; end
   hf = bar(time-dt,tJobkg,0.6); hold on
   ha = bar(time+dt,tJoana,0.6); hold off
   set(hf,'FaceColor','b','EdgeColor','none')
   set(ha,'FaceColor','r','EdgeColor','none')
   ylim = [0 1.1*max([tJobkg' tJoana'])];
   ptaxes(time,ylim,1,'Normalized cost:','k','Jo(O-B)/p','b','Jo(O-A)/p','r')

%------------------------------------------------------------------------------

function ptaxes(time,ylim,ticklabels,varargin)

date = floor(time);
hour = round((time-date)*24);

xlim = [time(1) time(end)] + 0.5*(time(end)-time(1))*[-1 1]/length(time);
set(gca,'XLim',xlim); 
set(gca,'YAxisLocation','right'); 
if diff(ylim), set(gca,'YLim',ylim); end 

set(gca,'XTick',time(hour==0))   % x-tick marks at 0Z
set(gca,'XTickLabel',[])         % x-tick mark labels: see below

if ticklabels,

   ith = find(hour==0);
   while length(ith)>10*2,          % thin until no more than 10 labels
      ith=ith(1:2:end);
   end

   dy = 0.02*diff(ylim);
%    text(time(ith(min(2,end))),ylim(2)+dy,'0Z',...
%                         'FontSize',8,...
%                         'HorizontalAlignment','center',...
%                         'VerticalAlignment','bottom');
   for i = ith,
     str = [datestr(date(i),'dd') datestr(date(i),'mmm')];
     if (date(i)==max(date)|i==ith(end)),    % last one: add year
        str = [str datestr(date(i),'yyyy')];
     end
     text(time(i),-dy,str,...
                        'FontSize',8,...
                        'HorizontalAlignment','left',...
                        'VerticalAlignment','top');
   end

end

% write title/legends:

set(gca,'YLimMode','manual')
if nargin>3,    
    np = round((nargin-3)/2);
    x = xlim(1); 
    y = ylim(2)+0.02*diff(ylim);
    for i = 1:np
        h = text(x,y,[' ' varargin{2*i-1}],'Color',varargin{2*i},...
                        'FontSize',8,'FontWeight','b',...
                        'HorizontalAlignment','left',...
                        'VerticalAlignment','bottom');
        pos = get(h,'Extent');
        x = pos(1) + pos(3);
    end   
end

%  add grid lines:

set(gca,'YGrid','on');
set(gca,'XGrid','on');

