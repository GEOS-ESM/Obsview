function ptotals(hfig,ts,ir,expid)

% PTOTALS Plot time-averaged statistics for OBS2HTML.

%     ts.id
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
%     ts.somfa(reg,lev,time)
%     ts.Jobkg(reg,lev,time)
%     ts.Joana(reg,lev,time)

if nargin<3, expid = ''; end

% sum over time:

nobs  = sum(ts.nobs (ir,:,:),3)';
nana  = sum(ts.nana (ir,:,:),3)';
npsv  = sum(ts.npsv (ir,:,:),3)';
somf  = sum(ts.somf (ir,:,:),3)';
somf2 = sum(ts.somf2(ir,:,:),3)';
soma  = sum(ts.soma (ir,:,:),3)';
soma2 = sum(ts.soma2(ir,:,:),3)';
somfa = sum(ts.somfa(ir,:,:),3)';
sigo  = sum(ts.sigo (ir,:,:),3)';
sigo2 = sum(ts.sigo2(ir,:,:),3)';
Jobkg = sum(ts.Jobkg(ir,:,:),3)';
Joana = sum(ts.Joana(ir,:,:),3)';

nexcl = nobs - (nana + npsv);

% create bar plots:

set(hfig,'Units','points','DefaultAxesFontSize',10)
figure(hfig), clf

if any(ts.lev),
  lev = ts.lev;
else
  lev = 1;
end
nlev = length(lev);
zero = zeros(nlev,1);
nans = NaN(nlev,1);

% Plot title
% ----------
subplot('Position',[0.1 0.9 0.8 0.1])

   t1 = ts.time(1); t2 = ts.time(end);
   ptitle = [datestr(t1,'dd') datestr(t1,'mmm') datestr(t1,'yyyy') ' ' datestr(t1,'HH') 'Z'];
   if t2>t1,
       ptitle = [ptitle ' - ' datestr(t2,'dd') datestr(t2,'mmm') datestr(t2,'yyyy') ' ' datestr(t2,'HH') 'Z'];
   end
   if any(expid), ptitle = [expid '  ' ptitle]; end
   text(0.5,0.5,ptitle,'FontSize',12,'FontWeight','b','Color','k',...
                       'HorizontalAlignment','center',...
                       'VerticalAlignment','bottom',...
                       'Interpreter','none');
   ptitle = [ts.id '  (' ts.reg{ir} ')'];
   text(0.5,0.1,ptitle,'FontSize',12,'FontWeight','b','Color','b',...
                       'HorizontalAlignment','center',...
                       'VerticalAlignment','bottom',...
                       'Interpreter','none');
                
   set(gca,'Visible','off')

% Data counts
% -----------
ph = min(0.35,nlev*0.1);
subplot('Position',[0.1 0.87-ph 0.35 ph])

   if nlev>1
       h = barh(1:nlev, [nana npsv nexcl],'stack');
   else
       h = barh([1 1 1; 0 0 0], [nana npsv nexcl; 0 0 0],'stack');
   end
   set(h,{'FaceColor'},{'g';[1 0.5 0];'r'})
   xlim = [0 1.1*max(nobs)];
   ptaxis(lev,xlim,'Data counts:','k','used (p)','g','+','k','passive',[1 0.5 0],'+','k','not used','r')

% Omf, oma statistics
% -------------------
subplot('Position',[0.6 0.87-ph 0.35 ph])

   n = nana + npsv;
   i = (n>0);
   meanf = zero; meanf(i) = somf(i)./n(i);
   meana = zero; meana(i) = soma(i)./n(i);
   rmsf  = zero; rmsf(i) = sqrt(somf2(i)./n(i));
   rmsa  = zero; rmsa(i) = sqrt(soma2(i)./n(i));

   hrf = barh((1:nlev)+0.1,rmsf,0.6); hold on
   hmf = barh((1:nlev)+0.1,meanf,0.6); hold on
   hra = barh((1:nlev)-0.1,rmsa,0.6); hold on
   hma = barh((1:nlev)-0.1,meana,0.6); hold off
   set(hrf,'FaceColor','b','EdgeColor','none')
   set(hmf,'FaceColor','c','EdgeColor','none')
   set(hra,'FaceColor','r','EdgeColor','none')
   set(hma,'FaceColor',[1 0.5 0],'EdgeColor','none')   
   xmax = max([rmsf' rmsa']);
   xmin = min([0 meanf' meana']);
   if xmin<0, xmin = min(xmin,-0.15*xmax); end
   xlim = 1.1*[xmin xmax];
   ptaxis(lev,xlim,'mean(O-B)','c','rms(O-B)','b','mean(O-A)',[1 0.5 0],'rms(O-A)','r')
    
% Jo(bkg)/n, Jo(ana)/n
% --------------------
subplot('Position',[0.1 0.75-2*ph 0.35 ph])

   n = nana;
   i = (n>0);
   Jobkg(i) = Jobkg(i)./n(i);
   Joana(i) = Joana(i)./n(i);
   hf = barh((1:nlev)+0.1,Jobkg,0.6); hold on
   ha = barh((1:nlev)-0.1,Joana,0.6); hold off
   set(hf,'FaceColor','b','EdgeColor','none')
   set(ha,'FaceColor','r','EdgeColor','none') 
   xlim = [0 1.1*max([Jobkg' Joana'])];
   ptaxis(lev,xlim,'Normalized cost:','k','Jo(O-B)/p','b','Jo(O-A)/p','r')

% Mean and rms sigo
% ------------------
subplot('Position',[0.6 0.75-2*ph 0.35 ph])

   n = nana;
   i = (n>0);
   tmean = zero; tmean(i) = sigo(i)./n(i);
   trms  = zero; trms(i)  = sqrt(sigo2(i)./n(i));
   
   n = nana + npsv;
   i = (n>0);
   esigo = nans; esigo(i) = sqrt(max(0,somfa(i))./n(i));
   esigb = nans; esigb(i) = sqrt(max(0,somf2(i)-somfa(i))./n(i));
   
   hr = barh(1:nlev,trms ); hold on
   hm = barh(1:nlev,tmean); hold on
   set(hr,'FaceColor','b')
   set(hm,'FaceColor','c')
   plot(esigo,(1:nlev)+0.1,'LineStyle','none',...
       'Marker','o','MarkerFaceColor',[1 0.5 0],'MarkerEdgeColor','k'); hold on
   plot(esigb,(1:nlev)-0.1,'LineStyle','none',...
       'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k'); hold off
   xlim = [0 1.1*max([trms' esigo' esigb'])];
   ptaxis(lev,xlim,'mean(sigo)','c','rms(sigo)','b','est(sigo)',[1 0.5 0],'est(sigb)','k')

%------------------------------------------------------------------------------

function ptaxis(lev,xlim,varargin) 

nlev = length(lev);
set(gca,'YLim',[0.5 nlev+0.5],...
        'YTick',[1:nlev],...
        'YTickLabel',lev)       
if diff(xlim), set(gca,'XLim',xlim); end
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');

set(gca,'XLimMode','manual')
np = round((nargin-2)/2);
x = xlim(1);
y = ylim(2)+0.02*diff(ylim);
for i = 1:np
    h = text(x,y,[' ' varargin{2*i-1}],'Color',varargin{2*i},...
        'FontSize',10,'FontWeight','b',...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom');
    pos = get(h,'Extent');
    x = pos(1) + pos(3);
end
