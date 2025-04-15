function ods = getods(varargin)

% GETODS Read raw data from an ODS file (obsolete: use ODSLOAD).
%
%    ODS = GETODS(FNAME) reads all data from the ODS file FNAME, including
%    file header information, and stores them in the structure ODS.
%
%    ODS = GETODS(FNAME,JDAY) returns data for Julian days in JDAY only.
%    ODS = GETODS(FNAME,JDAY,HOUR) returns data for Julian days in JDAY
%    and synoptic hours in HOUR only.  Note that JDAY and HOUR may be
%    scalars or vectors, but they must either be of the same length or
%    one of them must be a scalar.
%
%    Examples:
%
%    [fjday,ljday,lhour,nhour] = getodstimeinfo('file.ods');
%    ods = getods('file.ods',[fjday:ljday],0)
%    ods = getods('file.ods',fjday,[0 12])
%    ods = getods('file.ods',[fjday fjday+1],[0 12])
%
%    ODS = GETODS(FNAME,...,ATTR) will only return the ODS attributes
%    specified in the cell array ATTR.
%
%    Each ODS attribute will be stored as a column vector in a field
%    of the output structure. Each field's datatype will be inherited
%    from the file. If an offset and/or scale and/or missing value is
%    defined for an attribute, then this information will be added to
%    the output structure as well.
%
%    Use CONVERTODS to convert the output of GETODS to a structure
%    with double precision fields, properly translated and scaled,
%    and with missing values replaced by NaNs.
%
%    Use GETODSINFO to get file header information only.
%    Use GETODSTIMEINFO to get time information only.
%
%    This routine cannot read ODS prior to version 2.00
%
% This function is obsolete: Use ODSLOAD instead of GETODS.

warning('This function is obsolete: Use ODSLOAD instead of GETODS.')

% 25Mar2002 Dick Dee (dee@dao.gsfc.nasa.gov)

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

% first check the last argument (if any):

if nargs>0 && iscell(args{nargs}),   % selected attributes only
   attrs = args{nargs};
   nargs = nargs-1;
   [fjday,ljday,lhour,nhour] = getodstimeinfo(odsfile);
   ods = [];
else                                % all attributes
   attrs = {'kt','kx','ks','lon','lat','lev','time',...
            'obs','omf','oma','xm','qcexcl','qchist','xvec'};
   ods = getodsinfo(odsfile);     % complete header info
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

nsonfile = 1 + nhour*(ljday - fjday + lhour/24);
is = [];

switch nargs

case 0,   % load entire file

   is = 1:nsonfile;

case 1,   % data for one or more entire days

   jday = args{1};
   is = repmat(nhour*(jday(:)-fjday),[1 nhour])';
   is = is(:) + repmat(1:nhour,[1 length(jday)])';
   is(is>nsonfile) = [];

case 2,   % data for one or more synoptic hours

   jday = args{1}; hour = args{2};
   if length(jday)>1&&length(hour)>1&&length(jday)~=length(hour),
      error('Invalid time specification.')
   end
   is = 1 + nhour*((jday(:)-fjday) + hour(:)/24);
   is(is>nsonfile) = [];

end

if isempty(is), return, end

if any(is<1|is>nsonfile)||~isequal(is,round(is)),
   error('Invalid time specification.')
end

% combine synoptic hours where possible:

i = find(diff(is)~=1)';
isbeg = is(1+[0 i]);
isend = is([i length(is)]);

% get the data from file:

ods = getdata(ods,odsfile,isbeg,isend,attrs);

%--------------------------------------------------------------------

function ods = getdata(ods,odsfile,isbeg,isend,attrs)

% this function reads ODS attribute data from file

% open the file:

SD_id = hdfsd('start',odsfile,'read');
if SD_id==-1, error(['Can''t open ' odsfile]); end

% get pointers to synoptic times on file:

[idx,nob] = getstidx(SD_id);

% idx: column vector of indices of first obs for each synoptic time
% nob: column vector with total number of obs for each synoptic time

for k = 1:length(isbeg),

    ibeg = idx(isbeg(k));                           % index of first obs to be read
    iend = ibeg + sum(nob(isbeg(k):isend(k))) - 1;  % index of last obs to be read

    for name = attrs,

        dset_name = name{1};
        dset_idx = hdfsd('nametoindex',SD_id,dset_name);
        sds_id = hdfsd('select',SD_id, dset_idx);
        [dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);

        if rank==2,

           blen = dimsizes(2);                % length of each batch
           bbeg = floor((ibeg-1)/blen);       % index (zero-base) of first batch to read
           bend = floor((iend-1)/blen);       % index (zero-base) of last batch to read

           start  = [bbeg 0];
           stride = [1 1];
           edge   = [bend-bbeg+1 blen];
           [data,status] = hdfsd('readdata',sds_id,start,stride,edge);
           if status==-1,
              hdfsd('end',SD_id);
              error(['Can''t read ' dset_name '.'])
           end
           i1 = rem(ibeg-1,blen)+1; i2 = i1 + (iend-ibeg);
           attr = data(i1:i2);
           ods.(dset_name) = [ods.(dset_name); attr(:)];

           for attr_idx = 0:nattrs-1,

               [name,data_type,count,status] = hdfsd('attrinfo',sds_id,attr_idx);
               if status==-1, error('hdfsd(''attrinfo'') failed.'); end
               if strcmp(name,'add_offset'),
                  [add_offset,status] = hdfsd('readattr',sds_id,attr_idx);
                  if status==-1, error('hdfsd(''readattr'') failed.'); end
                  ods.([dset_name '_offset']) = add_offset;
               end
               if strcmp(name,'scale_factor'),
                  [scale_factor,status] = hdfsd('readattr',sds_id,attr_idx);
               if status==-1, error('hdfsd(''readattr'') failed.'); end
               ods.([dset_name '_scale']) = scale_factor;
               end
               if strcmp(name,'missing_value'),
                  [missing_value,status] = hdfsd('readattr',sds_id,attr_idx);
                  if status==-1, error('hdfsd(''readattr'') failed.'); end
                  ods.([dset_name '_missing']) = missing_value;
               end

           end

        end

        hdfsd('endaccess',sds_id);

    end

end

% close the file:

status = hdfsd('end',SD_id);
if status==-1, error(['Can''t close ' odsfile]); end

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
   idx = hdfsd('nametoindex',SD_id,dset_name);
   sds_id = hdfsd('select',SD_id, idx);
   [dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
   if status==-1,
        hdfsd('end',SD_id);
        error(['Can''t read ' dset_name '.'])
   end

   if all(isfinite(dimsizes)) && length(dimsizes)==2,

      start  = [0 0];
      stride = [1 1];
      edge   = dimsizes;
      [data,status] = hdfsd('readdata',sds_id,start,stride,edge);
      if status==-1,
         hdfsd('end',SD_id);
         error(['Can''t read ' dset_name '.'])
      end

      collect = [collect double(data(:))];
      hdfsd('endaccess',sds_id);

   else

      hdfsd('end',SD_id);
      error(['Unexpected dimensions of ' dset_name '.']);

   end

end

idx  = collect(:,1); % beginning indices for each synoptic time
nobs = collect(:,2); % number of obs for each synoptic time

n = sum(nobs);       % total number of obs on file
i = find(idx<=n);    % remove out-of-bound indices:
idx  = idx(i);
nobs = nobs(i);
