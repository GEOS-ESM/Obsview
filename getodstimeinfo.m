function [fjday,ljday,lhour,nsyn,varargout] = getodstimeinfo(odsfile)

% GETODSTIMEINFO Read time information from an ODS file.
%
%    [FJDAY,LJDAY,LHOUR,NSYN] = GETODSTIMEINFO(FNAME) returns the first Julian day
%    FJDAY, last Julian day LJDAY, last synoptic hour LHOUR, number of synoptic
%    hours per day NSYN, defined in the header of the ODS file FNAME.
%
%    [FJDAY,LJDAY,LHOUR,NSYN,NOBS] = GETODSTIMEINFO(FNAME) also returns the
%    array NOBS containing the number of observations for each of the synoptic
%    hours on file.

% 25Mar2002 Dick Dee (dee@dao.gsfc.nasa.gov)
% 18Jan2005 - added NOBS
% 15Apr2025 - updated to matlab 2024a
if nargin==0, help(mfilename), return, end

% first argument must be the name of a file:

if ~exist(odsfile,'file'),
    help(mfilename); error([odsfile ': No such file.']);
end

% open the file:

%SD_id = hdfsd('start',odsfile,'read');
SD_id = netcdf.open(odsfile);
%if SD_id==-1, error(['Can''t open ' odsfile]); end
% get time information:
%idx = hdfsd('nametoindex',SD_id,'syn_beg');
idx = netcdf.inqVarID(SD_id,'syn_beg');
%sds_id = hdfsd('select',SD_id, idx);
%[dname,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id);
[dname,data_type,dimsizes,nattrs] = netcdf.inqVar(SD_id,idx);
%if status==-1, warning('Can''t read ''syn_beg''.'); end
nsyn = dimsizes(2);

%attr_index = hdfsd('findattr',sds_id,'first_julian_day');
%[data,status] = hdfsd('readattr',sds_id,attr_index);
data = netcdf.getAtt(SD_id,idx,'first_julian_day');
%if status==-1, warning(['Can''t read ' attr_name '.']); end
fjday = double(data);

%attr_index = hdfsd('findattr',sds_id,'latest_julian_day');
data = netcdf.getAtt(SD_id,idx,'latest_julian_day');
%[data,status] = hdfsd('readattr',sds_id,attr_index);

%if status==-1, warning(['Can''t read ' attr_name '.']); end
ljday = double(data);

data = netcdf.getAtt(SD_id,idx,'latest_synoptic_hour');
%attr_index = hdfsd('findattr',sds_id,'latest_synoptic_hour');
%[data,status] = hdfsd('readattr',sds_id,attr_index);

%if status==-1, warning(['Can''t read ' attr_name '.']); end
lhour = double(data);
if nargout==5,
    disp('entering getstidx')
    [idx,nob] = getstidx(SD_id);
    disp('exiting getstidx')
    varargout(1) = {nob}

end

%hdfsd('endaccess',sds_id);

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
   [dimname2, dimlen2] = netcdf.inqDim(SD_id,dimids(2));
   [dimname1, dimlen1] = netcdf.inqDim(SD_id,dimids(1));
   dimsizes = [dimlen2 dimlen1];

   % it's not dimsizes it's dimids
   % need to pull dimids to get dimsizes
%   if status==-1,
%        hdfsd('end',SD_id);
%        error(['Can''t read ' dset_name '.'])
%   end

   if all(isfinite(dimsizes)) && length(dimsizes)==2,

      start  = [0 0];
      stride = [1 1];
      edge   = dimsizes;
      dname
      data = netcdf.getVar(SD_id,idx); %,start,stride,edge)
      %[data,status] = hdfsd('readdata',sds_id,start,stride,edge);
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
n = sum(nobs);       % total number of obs on file
idx;
i = find(idx<=n);    % remove out-of-bound indices:
idx  = idx(i);
%nobs = sum(nobs(i))
nobs = nobs(i)
