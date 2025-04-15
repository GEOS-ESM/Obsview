function ods = setcidx(ods,coloring) 

% SETCIDX Set color index attribute in ODS structure.

ods.cidx = ones(size(ods.kt),'int8');
ods.cinfo = [];

switch coloring.cmap

    case 'bwr'         % create a blue-white-red colormap

        a = linspace(0,1,32)';
        b = linspace(1,0,32)';
        c = ones(32,1);
        cmap = [a a c; c b b];

    otherwise

        cmap = colormap(coloring.cmap);

end

nr = size(cmap,1);
scol = [0.8125 1.0000 0.1875];
scol = [0 0 1];

attr = coloring.attr;

switch attr

   case 'kt'

      kts = sort(unique(ods.kt));
      nc = length(kts);
      if nc==1,
         ods.cidx = ones(size(ods.kt),'int8');
         ods.cinfo(1).rgb = scol;
         ods.cinfo(1).txt = int2str(kts);
      else
         for i = 1:nc,
            ods.cidx(ods.kt==kts(i)) = i;
            ic = round(1 + (i-1)*(nr-1)/(nc-1));
            ods.cinfo(i).rgb = cmap(ic,:);
            ods.cinfo(i).txt = int2str(kts(i));
         end
      end

   case 'kx'

      kxs = sort(unique(ods.kx));
      nc = length(kxs);
      if nc==1,
         ods.cidx = ones(size(ods.kt),'int8');
         ods.cinfo(1).rgb = scol;
         ods.cinfo(1).txt = int2str(kxs);
      else
         for i = 1:nc,
            ods.cidx(ods.kx==kxs(i)) = i;
            ic = round(1 + (i-1)*(nr-1)/(nc-1));
            ods.cinfo(i).rgb = cmap(ic,:);
            ods.cinfo(i).txt = ['kx = ' int2str(kxs(i))];
         end
      end

   case 'lat'

      lats = -90:15:90;
      nc = length(lats)-1;
      for i = 1:nc-1,
         ods.cidx(lats(i)<=ods.lat&ods.lat<lats(i+1)) = i;
         ic = round(1 + (i-1)*(nr-1)/(nc-1));
         ods.cinfo(i).rgb = cmap(ic,:);
         ods.cinfo(i).txt = ['[' latstr(lats(i)) ',' latstr(lats(i+1)) ')'];
      end
      i = nc;
      ods.cidx(lats(i)<=ods.lat&ods.lat<=lats(i+1)) = i;
      ic = round(1 + (i-1)*(nr-1)/(nc-1));
      ods.cinfo(i).rgb = cmap(ic,:);
      ods.cinfo(i).txt = ['[' latstr(lats(i)) ',' latstr(lats(i+1)) ']'];

   case 'lon'

      lons = -180:30:180;
      nc = length(lons)-1;
      for i = 1:nc,
         ods.cidx(lons(i)<=ods.lon&ods.lon<lons(i+1)) = i;
         ic = round(1 + (i-1)*(nr-1)/(nc-1));
         ods.cinfo(i).rgb = cmap(ic,:);
         ods.cinfo(i).txt = ['[' lonstr(lons(i)) ',' lonstr(lons(i+1)) ')'];
      end

   case 'lev'

      ncmax = 20;
      p = sort(unique(ods.lev),'descend');
      nc = length(p);
      if nc==1,
         ods.cidx = ones(size(ods.kt),'int8');
         ods.cinfo(1).rgb = scol;
         ods.cinfo(1).txt = [num2str(p) 'hPa'];
      elseif nc<ncmax,  % discrete set of levels
         for i = 1:nc,
            ods.cidx(ods.lev==p(i)) = i;
            ic = round(1 + (i-1)*(nr-1)/(nc-1));
            ods.cinfo(i).rgb = cmap(ic,:);
            ods.cinfo(i).txt = [num2str(p(i)) 'hPa'];
         end
      else               % uniformly spaced bins (in log(p))
         nc = ncmax;
         logp = log(ods.lev);
         maxlogp = max(logp); minlogp = min(logp);
         dz = (maxlogp-minlogp)/nc;
         zz = minlogp + dz*(0:nc);
         zz(end) = maxlogp;
         zbins = zz + max(eps,eps*abs(zz));
         pbins = exp(zbins);
         pbins(1) = 0; pbins(end) = Inf;     % bin edges [ )
         for i = 1:nc,
            ods.cidx(pbins(nc-i+1)<=ods.lev&ods.lev<pbins(nc-i+2)) = i;
            ic = round(1 + (i-1)*(nr-1)/(nc-1));
            ods.cinfo(i).rgb = cmap(ic,:);
            if i==nc,
               ods.cinfo(i).txt = ['< ' num2str(pbins(nc-i+2),3) 'hPa'];
            elseif i==1,
               ods.cinfo(i).txt = ['> ' num2str(pbins(nc-i+1),3) 'hPa'];
            else
               ods.cinfo(i).txt = ['[' num2str(pbins(nc-i+1),3) 'hPa,' num2str(pbins(nc-i+2),3) 'hPa)'];
            end
         end
      end

   case {'time','obs','sigo','omf','xm','ks'}

      nc = 20;
      maxt = double(max(ods.(attr)));
      mint = double(min(ods.(attr)));
      dz = (maxt-mint)/nc;
      zz = mint + dz*(0:nc);
      zz(end) = maxt;
      zbins = zz + max(eps,eps*abs(zz));
      for i = 1:nc,
         ods.cidx(zbins(i)<=ods.(attr) & ods.(attr)<zbins(i+1)) = i;
         ic = round(1 + (i-1)*(nr-1)/(nc-1));
         ods.cinfo(i).rgb = cmap(ic,:);
         ods.cinfo(i).txt = [attr ' in [ ' num2str(zbins(i)) ', ' num2str(zbins(i+1)) ')'];
      end

   case 'qcx'

      ods.cidx(ods.qcx==0) = 1;
      ods.cinfo(1).rgb = [0 1 0];   % green
      ods.cinfo(1).txt = 'Used';
      ods.cidx(ods.qcx==1) = 2;
      ods.cinfo(2).rgb = [1 0.5 0];   % orange
      ods.cinfo(2).txt = 'Passive';
      ods.cidx(ods.qcx> 1) = 3;
      ods.cinfo(3).rgb = [1 0 0];   % red
      ods.cinfo(3).txt = 'Rejected by QC';

   case 'qch'

      qchs = sort(unique(ods.qch));
      nc = length(qchs);
      if nc==1,
         ods.cidx = ones(size(ods.qch),'int8');
         ods.cinfo(1).rgb = scol;
         ods.cinfo(1).txt = ['qch = ' int2str(qchs)];
      else
         for i = 1:nc,
            ods.cidx(ods.qch==qchs(i)) = i;
            ic = round(1 + (i-1)*(nr-1)/(nc-1));
            ods.cinfo(i).rgb = cmap(ic,:);
            ods.cinfo(i).txt = ['qch = ' int2str(qchs(i))];
         end
      end

   otherwise

      ods.cidx = ones(size(ods.kt),'int8');
      ods.cinfo.rgb = scol;
      ods.cinfo.txt = '';

end
