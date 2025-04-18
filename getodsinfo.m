function ods = getodsinfo(odsfile)
% GETODSINFO Read header information from an ODS file.
%
%    ODS = GETODSINFO(FNAME) returns all header information
%    from the ODS file FNAME.

% 25Mar2002 Dick Dee (dee@dao.gsfc.nasa.gov)
% 15Apr2025 Wesley Davis (wesley.j.davis@nasa.gov)
if nargin==0, help(mfilename), return, end

% first argument must be the name of a file:

if ~exist(odsfile,'file'),
    help(mfilename); error([odsfile ': No such file.']);
end

% try to open the file:
SD_id = netcdf.open(odsfile);

% return the file name in the output structure:

ods.filename = odsfile;

% get header information and return in the output structure:
ods = getinfo(ods,SD_id);

% close the file:

netcdf.close(SD_id);

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

varids=netcdf.inqVarIDs(SD_id);
ndatasets=length(varids);
finfo=ncinfo(ods.filename);
nglobal_attr=length(finfo.Attributes);
attrNames = {finfo.Attributes.Name};

% get all global attributes:

for attr_idx = 1:nglobal_attr,
    glbl_attr_name_cell  = attrNames(attr_idx);
    glbl_attr_name       = glbl_attr_name_cell{1};
    data                 = ncreadatt(ods.filename,'/',glbl_attr_name);
    if isnumeric(data), data = double(data); end
    ods.(glbl_attr_name) = data;

end

% get time information for data on file:
% pull only synoptic begin

idx = netcdf.inqVarID(SD_id,'syn_beg');
Varname=netcdf.inqVar(SD_id,idx);
sds_id = netcdf.getVar(SD_id,idx);

for name = {'first_julian_day','latest_julian_day','latest_synoptic_hour'},

    attr_name = char(name);
    data = netcdf.getAtt(SD_id,idx,attr_name);

    if isnumeric(data),
       data = double(data);
    else
       warning(['Unexpected numeric attribute: ' attr_name])
    end
    ods.(attr_name) = data;

end

[dname,xtype,dimids,nattrs] = netcdf.inqVar(SD_id,idx);
[dimname2, dimlen2] = netcdf.inqDim(SD_id,dimids(2));
[dimname1, dimlen1] = netcdf.inqDim(SD_id,dimids(1));
ods.synoptic_hours_per_day = dimlen1;

% get additional header information:

for name = {'kt_names','kt_units','kx_names','kx_meta','qcx_names'},

    dset_name = char(name);
    idx = netcdf.inqVarID(SD_id,dset_name);
    [dname,data_type,dimids,nattrs] = netcdf.inqVar(SD_id,idx);
    [dimname2, dimlen2] = netcdf.inqDim(SD_id,dimids(2));
    [dimname1, dimlen1] = netcdf.inqDim(SD_id,dimids(1));
    dimsizes = [dimlen2 dimlen1];

    if all(isfinite(dimsizes)),

       data = netcdf.getVar(SD_id,idx);
       ods.(dset_name) = data;

    end
    
end
