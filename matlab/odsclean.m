function ods = odsclean(ods)

% ODSCLEAN Clean up ODS structure: empty attributes, NaNs etc.

[MAXR,KTS] = dconfig('MAXR','KTS');

% replace empty attributes by NaNs:
% --------------------------------
nans = NaN*ones(size(ods.kt));
fields = fieldnames(ods)
for j = 1:length(fields),
    if isempty(ods.(fields{j})), ods.(fields{j}) = nans; end
end

% For obs2html to work in the absence of all oma information:
% ----------------------------------------------------------
if all(isnan(ods.oma))
   ods.oma = zeros(size(ods.oma)); 
end

% set exclusion flag for data with very large sigo's:
% --------------------------------------------------
for j = 1:length(KTS)
   i = (ods.kt==KTS(j).value & ods.sigo>KTS(j).msigo & ods.qcx==0); 
   ods.qcx(i) = 3;
end

% replace numerical undef's by NaNs:
% ---------------------------------
i = abs(ods.lev )>MAXR; ods.lev (i) = NaN;
i = abs(ods.obs )>MAXR; ods.obs (i) = NaN;
i = abs(ods.omf )>MAXR; ods.omf (i) = NaN;
i = abs(ods.oma )>MAXR; ods.oma (i) = NaN;
i = abs(ods.xm  )>MAXR; ods.xm  (i) = NaN;
i = abs(ods.sigo)>MAXR; ods.sigo(i) = NaN;

% set exclusion flag for data with undefined omf/oma:
% --------------------------------------------------
i = (isnan(ods.omf)|isnan(ods.oma) & (ods.qcx==0|ods.qcx==1)); 
ods.qcx(i) = 4;
