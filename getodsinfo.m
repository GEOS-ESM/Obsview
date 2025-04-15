function ods = getodsinfo(odsfile)

% GETODSINFO Read header information from an ODS file.
%
%    ODS = GETODSINFO(FNAME) returns all header information
%    from the ODS file FNAME.

% 25Mar2002 Dick Dee (dee@dao.gsfc.nasa.gov)

if nargin==0, help(mfilename), return, end

% first argument must be the name of a file:

if ~exist(odsfile,'file'),
    help(mfilename); error([odsfile ': No such file.']);
end

% try to open the file:

SD_id = hdfsd('start',odsfile,'read');
if SD_id==-1,
   status = hdfsd('end',SD_id);
   error(['Can''t open ' odsfile]);
end

% return the file name in the output structure:

ods.filename = odsfile;

% get header information and return in the output structure:

ods = getinfo(ods,SD_id);

% close the file:

status = hdfsd('end',SD_id);
if status==-1, error(['Can''t close ' odsfile]); end

% all done.

%--------------------------------------------------------------------

function ods = getinfo(ods,SD_id)

% this function gets header information from the ODS file and adds
% it to the structure.

% header information includes:
% - all global file attributes, whatever they are
% - synoptic time information
% - various character arrays

% get file info:

[ndatasets,nglobal_attr,status] = hdfsd('fileinfo',SD_id);
if status==-1, error('hdfsd(''fileinfo'') failed.'); end

% get all global attributes:

for attr_idx = 0:nglobal_attr-1,

    [attr_name,data_type,count,status] = hdfsd('attrinfo',SD_id,attr_idx);
    [data,status] = hdfsd('readattr',SD_id,attr_idx);
    if status==-1,
       hdfsd('end',SD_id);
       error('Can''t read global attributes.')
    end
    if isnumeric(data), data = double(data); end

    ods.(attr_name) = data;

end

% get time information for data on file:

idx = hdfsd('nametoindex',SD_id,'syn_beg');
sds_id = hdfsd('select',SD_id, idx);
if sds_id==-1,
   hdfsd('end',SD_id);
   error('Not an ODS file (Can''t read ''syn_beg'').')
end

for name = {'first_julian_day','latest_julian_day','latest_synoptic_hour'},

    attr_name = char(name);
    attr_index = hdfsd('findattr',sds_id,attr_name);
    [data,status] = hdfsd('readattr',sds_id,attr_index);
    if status==-1, warning(['Can''t read ' attr_name '.']); end

    if isnumeric(data),
       data = double(data);
    else
       warning(['Unexpected numeric attribute: ' attr_name])
    end

    ods.(attr_name) = data;

end

[dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
if status==-1,
   hdfsd('end',SD_id);
   error('Can''t read ''syn_beg''.')
end
ods.synoptic_hours_per_day = dimsizes(2);

hdfsd('endaccess',sds_id);

% get additional header information:

for name = {'kt_names','kt_units','kx_names','kx_meta','qcx_names'},

    dset_name = char(name);
    idx = hdfsd('nametoindex',SD_id,dset_name);
    sds_id = hdfsd('select',SD_id, idx);
    [dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
    if status==-1, warning(['Can''t read ' dset_name '.']); end

    if all(isfinite(dimsizes)),

       start  = zeros(size(dimsizes))
       stride = ones(size(dimsizes))
       edge   = dimsizes
       disp(start)
       disp(stride)
       disp(edge)
       [data,status] = hdfsd('readdata',sds_id,start,stride,edge);

       ods.(dset_name) = data;

    end

    hdfsd('endaccess',sds_id);

end

