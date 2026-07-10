function ods = convertods(ods)

% CONVERTODS Convert output from GETODS (obsolete: use ODSLOAD).

warning('This function is obsolete: Use ODSLOAD instead of GETODS.')

for attr = {'kt','kx','ks','lon','lat','lev','time',...
            'obs','omf','oma','xm','qcexcl','qchist','xvec'}

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
