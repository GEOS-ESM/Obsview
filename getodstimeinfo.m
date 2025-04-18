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

SD_id = netcdf.open(odsfile);
idx = netcdf.inqVarID(SD_id,'syn_beg');
% [dname,data_type,dimids,nattrs] = netcdf.inqVar(SD_id,idx);
% [dimname2, dimlen2] = netcdf.inqDim(SD_id,dimids(2));
% [dimname1, dimlen1] = netcdf.inqDim(SD_id,dimids(1));
% dimsizes = [dimlen2 dimlen1];
[dname,data_type,dimsizes,nattrs] = netcdf.inqVar(SD_id,idx);
nsyn = dimsizes(2);

data = netcdf.getAtt(SD_id,idx,'first_julian_day');
fjday = double(data);

data = netcdf.getAtt(SD_id,idx,'latest_julian_day');
ljday = double(data);

data = netcdf.getAtt(SD_id,idx,'latest_synoptic_hour');
lhour = double(data);

if nargout==5,
    [idx,nob] = getstidx(SD_id);
    varargout(1) = {nob}
end

% close the file:
netcdf.close(SD_id);
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
   [dname,data_type,dimids,nattrs] = netcdf.inqVar(SD_id,idx);
   [dimname2, dimlen2] = netcdf.inqDim(SD_id,dimids(2));
   [dimname1, dimlen1] = netcdf.inqDim(SD_id,dimids(1));
   dimsizes = [dimlen2 dimlen1];

   if all(isfinite(dimsizes)) && length(dimsizes)==2,

      data = netcdf.getVar(SD_id,idx);
      collect = [collect double(data(:))];

   else
      netcdf.close(SD_id);

   end

end

idx  = collect(:,1); % beginning indices for each synoptic time
nobs = collect(:,2); % number of obs for each synoptic time
n = sum(nobs);       % total number of obs on file
idx;
i = find(idx<=n);    % remove out-of-bound indices:
idx  = idx(i);
nobs = nobs(i)
