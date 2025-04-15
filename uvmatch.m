function [iu,iv] = uvmatch(ods)

% UVMATCH Find matching u-v wind components.
%
% [IU,IV] = UVMATCH(ODS) returns index vectors IU, IV into the
%       ods structure ODS for matching wind components. This means
%       that the kx, time, lon, lat, lev attributes must agree. 

% 25Oct2001 Dick Dee (dee@dao.gsfc.nasa.gov)
% 07Dec2004 Dick Dee - GSI version

[KTWND,KTPRS] = dconfig('KTWND','KTPRS'); 
    
iu = []; iv = [];

for ip = 1:length(KTWND),
    
    ktu   = KTWND{ip}(1);
    ktv   = KTWND{ip}(2);
    isprs = ismember(ktu,KTPRS);
    iup   = find(ods.kt==ktu);    % indices of all u-components...
    done  = 0;
    
    % first try if matching (u,v) always occur together (fast):
    
    ivp = iup + 1;
    if max(ivp)<=length(ods.kt) & all(ods.kt(ivp)==ktv),  % u,v appear in pairs
        if isequal(ods.kx(ivp),ods.kx(iup)) & ...
           isequal(ods.time(ivp),ods.time(iup)) & ...
           isequal(ods.lon(ivp),ods.lon(iup)) & ...
           isequal(ods.lat(ivp),ods.lat(iup)) & ...
           (isprs & isequal(ods.lev(ivp),ods.lev(iup))), done = 1;
        end
    end
    
    if ~done, % unfortunately we have to do a full search (slow)
        ivp = NaN*ones(size(iup));
        for k = 1:length(iup),     % for each u-component, find index of v-component
            i = iup(k);
            j = find(ods.kt==ktv & ods.kx==ods.kx(i) & ods.time==ods.time(i));   % from same instrument & sounding
            if any(j)&isprs, j = j(ods.lev(j)==ods.lev(i)); end  % .. at same level
            if any(j), j = j(ods.lon(j)==ods.lon(i) & ods.lat(j)==ods.lat(i)); end  % ..and same location
            if length(j)==1, ivp(k) = j; end  % we should have exactly one, else undefined
        end
    end
    
    iu = [iu; iup(isfinite(ivp))];
    iv = [iv; ivp(isfinite(ivp))];
    
end