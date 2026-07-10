function ods = odsload(varargin)
% ODSLOAD Load ODS data structure from file.
%
%    ODS = ODSLOAD(FNAME) reads all data from the ODS file FNAME, including
%    file header information, and stores them in the structure ODS.
%
%    ODS = ODSLOAD(FNAME,JDAY) returns data for Julian days in JDAY only.
%
%    ODS = ODSLOAD(FNAME,JDAY,HOUR) returns data for Julian days in JDAY
%    and synoptic hours in HOUR only.  Note that JDAY and HOUR may be
%    scalars or vectors, but they must either be of the same length or
%    one of them must be a scalar.
%
%    ODS = ODSLOAD(FNAME,...,ATTR) will only return the ODS attributes
%    specified in the cell array ATTR.
%
%    Examples:
%
%    [fjday,ljday,lhour,nhour] = getodstimeinfo('file.ods');
%    ods = odsload('file.ods',[fjday:ljday],0)
%    ods = odsload('file.ods',fjday,[0 12])
%    ods = odsload('file.ods',[fjday fjday+1],[0 12])
%    ods = odsload('file.ods',{'lat','lon'})
%
%    Each ODS attribute will be stored as a column vector in a field
%    of the output structure. Each field's datatype will be inherited
%    from the file. If an offset and/or scale and/or missing value is
%    defined for an attribute, then this information will be added to
%    the output structure as well.
%
%    Use GETODSINFO to get file header information only.
%    Use GETODSTIMEINFO to get time information only.

% 01Jul2005 Dick Dee (dee@gmao.gsfc.nasa.gov)

if nargin==0,           % must have at least one argument
   help(mfilename)
   return,
end

odsfile = varargin{1};  % first argument must be ODS file

if nargin>1,            % additional arguments. if any
   nargs = nargin-1;
   args = varargin(2:nargin);
else
   nargs = 0;
end
nargs
args
% first check the last argument (if any):

if nargs>0 && iscell(args{nargs}),   % selected attributes only
   attrs = args{nargs};
   for i = 1:length(attrs) % map GSI attributes to PSAS attributes
      if strcmp(attrs{i},'sigo'), attrs{i} = 'xvec'  ; end
      if strcmp(attrs{i},'qch' ), attrs{i} = 'qchist'; end
      if strcmp(attrs{i},'qcx' ), attrs{i} = 'qcexcl'; end
   end
   nargs = nargs-1;
   [fjday,ljday,lhour,nhour] = getodstimeinfo(odsfile);
   ods = [];
else                                % all attributes
   attrs = {'kt','kx','ks','lon','lat','lev','time',...
            'obs','omf','oma','xm','qcexcl','qchist','xvec'};
   disp('enter getodsinfo from odsload')
   ods = getodsinfo(odsfile)     % complete header info
   disp('exiting getodsinfo from odsload')
   if ~(isodsstruct(ods)&&str2num(ods.version(1))>=2),
      error([odsfile ': Not an ODS Version 2 file.'])
   end
   fjday = ods.first_julian_day;
   ljday = ods.latest_julian_day;
   lhour = ods.latest_synoptic_hour;
   nhour = ods.synoptic_hours_per_day;
end
% initialize attribute fields:

for field = attrs,
   ods.(field{1}) = [];
end

% any remaining arguments specify synoptic times:

nsonfile = 1 + nhour*(ljday - fjday + lhour/24)
is = []
disp('switch nargs in odsload')
switch nargs

case 0,   % load entire file

   is = 1:nsonfile;
   disp('nargs = 0')

case 1,   % data for one or more entire days

   jday = args{1};
   is = repmat(nhour*(jday(:)-fjday),[1 nhour])';
   is = is(:) + repmat(1:nhour,[1 length(jday)])';
   is(is>nsonfile) = [];
   disp('nargs = 1')

case 2,   % data for one or more synoptic hours

   jday = args{1} 
   hour = args{2}
   if length(jday)>1&&length(hour)>1&&length(jday)~=length(hour),
      error('Invalid time specification.')
   end
   is = 1 + nhour*((jday(:)-fjday) + hour(:)/24)
   is(is>nsonfile) = []
   disp('nargs = 2')

end
is
if isempty(is), return, end

if any(is<1|is>nsonfile)||~isequal(is,round(is)),
   error('Invalid time specification.')
end

% combine synoptic hours where possible:

i = find(diff(is)~=1)';
isbeg = is(1+[0 i]);
isend = is([i length(is)]);

% get the raw data from file:
disp('entering getdata in odsload')
ods = getdata(ods,odsfile,isbeg,isend,attrs);
disp('exiting getdata in odsload')
% scale and offset:
disp('entering convert in odsload')
ods = convert(ods,attrs);
disp('entering convert in odsload')
ods
%--------------------------------------------------------------------

function ods = getdata(ods,odsfile,isbeg,isend,attrs)

% this function reads ODS attribute data from file
% open the file:
SD_id = netcdf.open(odsfile);
%SD_id = hdfsd('start',odsfile,'read');
%if SD_id==-1, error(['Can''t open ' odsfile]); end

% get pointers to synoptic times on file:
[idx,nob] = getstidx(SD_id);
% idx: column vector of indices of first obs for each synoptic time
% nob: column vector with total number of obs for each synoptic time

for k = 1:length(isbeg),
    ibeg = idx(isbeg(k));                          % index of first obs to be read
    iend = ibeg + sum(nob(isbeg(k):isend(k))) - 1;  % index of last obs to be read

    for name = attrs,

        dset_name = char(name{1});
	
	dset_idx = netcdf.inqVarID(SD_id,dset_name);
	[dname,data_type,dimids,nattrs] = netcdf.inqVar(SD_id,dset_idx);
	[dimname1,dimsize1] = netcdf.inqDim(SD_id,dimids(1));
        [dimname2,dimsize2] = netcdf.inqDim(SD_id,dimids(2));
	
        %dset_idx = hdfsd('nametoindex',SD_id,dset_name);

        %sds_id = hdfsd('select',SD_id, dset_idx);
        %[dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);

        if length(dimids)==2,

           [dimname blen] = netcdf.inqDim(SD_id,dimids(1));                % length of each batch
           bbeg = floor((ibeg-1)/blen);       % index (zero-base) of first batch to read
           bend = floor((iend-1)/blen);       % index (zero-base) of last batch to read

           start  = [bbeg 0];
           stride = [1 1];
           edge   = [bend-bbeg+1 blen];
	   data = netcdf.getVar(SD_id,dset_idx); %,start,stride,edge)
	   
           %[data,status] = hdfsd('readdata',sds_id,start,stride,edge);
           %if status==-1,
           %   hdfsd('end',SD_id);
           %   error(['Can''t read ' dset_name '.'])
           %end
           i1 = rem(ibeg-1,blen)+1;
	   i2 = i1 + (iend-ibeg);
           attr = data(i1:i2);
           ods.(dset_name) = [ods.(dset_name); attr(:)];
           for attr_idx = 0:nattrs-1,
	       name = netcdf.inqAttName(SD_id,dset_idx,attr_idx);

               %[name,data_type,count,status] = hdfsd('attrinfo',sds_id,attr_idx);
               %if status==-1, error('hdfsd(''attrinfo'') failed.'); end
               if strcmp(name,'add_offset'),
		  add_offset = netcdf.getAtt(SD_id,dset_idx,name);
                  %[add_offset,status] = hdfsd('readattr',sds_id,attr_idx);
                  %if status==-1, error('hdfsd(''readattr'') failed.'); end
                  ods.([dset_name '_offset']) = add_offset;
               end
               if strcmp(name,'scale_factor'),
                  %[scale_factor,status] = hdfsd('readattr',sds_id,attr_idx);
		  scale_factor = netcdf.getAtt(SD_id,dset_idx,name);
            	  %if status==-1, error('hdfsd(''readattr'') failed.'); end
                  ods.([dset_name '_scale']) = scale_factor;
               end
               if strcmp(name,'missing_value'),
		  missing_value = netcdf.getAtt(SD_id,dset_idx,name);
                  %[missing_value,status] = hdfsd('readattr',sds_id,attr_idx);
                  %if status==-1, error('hdfsd(''readattr'') failed.'); end
                  ods.([dset_name '_missing']) = missing_value;
               end

           end

        end

        %hdfsd('endaccess',sds_id);

    end
end
% close the file:
netcdf.close(SD_id);
%status = hdfsd('end',SD_id);
%if status==-1, error(['Can''t close ' odsfile]); end

%--------------------------------------------------------------------

function [idx,nobs] = getstidx(SD_id);
% this function returns:
%
%      idx : column vector of indices of first obs for each synoptic time
%      nobs: column vector with total number of obs for each synoptic time

% synoptic time indices are stored in datasets called 'syn_beg','syn_len'

collect = [];

for name = {'syn_beg','syn_len'},
   dset_name = char(name);
   idx = netcdf.inqVarID(SD_id,dset_name);
   %idx = hdfsd('nametoindex',SD_id,dset_name);
   %sds_id = hdfsd('select',SD_id, idx);
   %[dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
   [dname,data_type,dimids,nattrs] = netcdf.inqVar(SD_id,idx);
   % this actually returns dimension IDs, the old one returns dimension sizes
   % so search for each dimension ID and return the names to see what they are
   [dimname2, dimlen2] = netcdf.inqDim(SD_id,dimids(2));
   
   [dimname1, dimlen1] = netcdf.inqDim(SD_id,dimids(1));
   dimsizes = [dimlen2 dimlen1];
   %if status==-1,
   %     hdfsd('end',SD_id);
   %     error(['Can''t read ' dset_name '.'])
   %end
   %ncid = netcdf.open("example.nc","NOWRITE");
   
   % Retrieve identifier of dimension.
   %dimid = netcdf.inqDimID(ncid,dimname)
   if all(isfinite(dimsizes)) && length(dimids)==2,

      start  = [0 0];
      stride = [1 1];
      edge   = dimsizes;
      %[data,status] = hdfsd('readdata',sds_id,start,stride,edge);
      %data = netcdf.getVar(SD_id,idx,start,stride,edge)
      data = netcdf.getVar(SD_id,idx);
      %if status==-1,
      %   hdfsd('end',SD_id);
      %   error(['Can''t read ' dset_name '.'])
      %end

      collect = [collect double(data(:))];
      %hdfsd('endaccess',sds_id);

   else

      %hdfsd('end',SD_id);
      netcdf.close(SD_id);
      error(['Unexpected dimensions of ' dset_name '.']);

   end

end

idx  = collect(:,1); % beginning indices for each synoptic time
nobs = collect(:,2); % number of obs for each synoptic time
n = sum(nobs);% total number of obs on file
i = find(idx<=n);    % remove out-of-bound indices:
idx  = idx(i)
nobs = nobs(i)
%--------------------------------------------------------------------

function ods = convert(ods,attrs)
% applies scaling and offsets, and changes some attribute names

for attr = attrs

    field = attr{:};
    if isfield(ods,field),
        
        data = ods.(field);
        ods = rmfield(ods,field);
                   
        modify = 0;
        field_offset = [field '_offset'];
        if isfield(ods,field_offset),
            offset = double(ods.(field_offset));
            ods = rmfield(ods,field_offset);
            if offset~=0, modify = 1; end
        else, offset = 0; end
        field_scale = [field '_scale'];
        if isfield(ods,field_scale),
            scale = double(ods.(field_scale));
            ods = rmfield(ods,field_scale);
            if scale~=1, modify = 1; end
        else, scale = 1; end
        field_missing = [field '_missing'];
        if isfield(ods,field_missing),
            modify = 1;
            missing = double(ods.(field_missing));
            ods = rmfield(ods,field_missing);
        else, missing = NaN; end
        
        if modify, 
            data = double(data);
            data = offset + scale*data;
            if isfinite(missing), data(data==missing) = NaN; end
        end

        switch field
            case 'kt'    
                ods.(field) = int8(data);
            case 'qcexcl'
                ods.qcx = int8(data);
            case {'kx'}
                ods.(field) = int16(data);
            case 'qchist'
                ods.qch = int16(data);
            case 'ks'
                ods.(field) = int32(data);
            case {'lon','lat','lev','time','obs','omf','oma','xm'}
                ods.(field) = single(data);  
            case 'xvec'
                ods.sigo = single(data); 
        end
               
    end
    
end
