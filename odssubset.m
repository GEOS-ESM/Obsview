function ods = odssubset(ods,arg)

% ODSSUBSET Return a subset of an ODS structure.
%
%    ODS = ODSSUBSET(ODS,I) returns the subset of the input ODS
%          structure indicated by the index (or logical) vector i. 
%    ODS = ODSSUBSET(ODS,GROUP) returns the subset defined by 
%          the attributes of GROUP: see OBSGROUPS. 

% 28Mar01 Dick Dee (dee@dao.gsfc.nasa.gov)
% 10Jan05 Generalized

attr = dconfig('OBSATTRIBUTES');
nsiz = size(ods.kt);

if isstruct(arg),

    grp = arg;
    is = true(nsiz);

    for j = 1:length(attr.names)
        a = attr.names{j};
        if isfield(grp,a)&&~isempty(grp.(a))  % select on this attribute
            if attr.discr(j)  % discrete case (list)
                i = false(size(is));
                for val = grp.(a), i = i|ods.(a)==val; end
                is = is & i;
            else         % continuous case (range, or single value)
                if strcmp(a,'lon'),  % special case (periodic)
                    lon = -180 + mod(grp.lon + 180,360);   % lon in [-180,180)
                    if length(lon)>1,
                        if lon(1)<lon(2),
                            is = is & (lon(1)<=ods.lon & ods.lon<lon(2));
                        elseif lon(1)>lon(2),
                            is = is & (lon(1)<=ods.lon | ods.lon<lon(2));
                        end
                    else
                        is = is & ods.lon==lon;
                    end
                else
                    if length(grp.(a))>1,
                        is = is & min(grp.(a))<=ods.(a) & ods.(a) < max(grp.(a));
                    else
                        is = is & ods.(a)==grp.(a);
                    end
                end
            end
        end
    end
    
else

    is = arg;

end

% subset all compatible fields

names = fieldnames(ods);
for i = 1:length(names)
    a = names{i};
    if isequal(size(ods.(a)),nsiz),
        ods.(a) = ods.(a)(is);
    end
end

