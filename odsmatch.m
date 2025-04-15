function x = odsmatch(x,y,fm,fc)

% ODSMATCH Match and copy fields of two data structures.
%
% For example, 
%
% ods = odsmatch(ods,ods2,{'kt','kx','lat','lon','lev','time'},{'obs','omf'})
%
% Just try it!

nm = length(fm); % fields to match
nc = length(fc); % fields to copy when a match is found

% sort x and y by fields to match

nx = numel(x.(fm{1}));
is = 1:nx;
for k = 1:nm
    f = fm{k};
    if numel(x.(f))~=nx, error('Wrong field length.'), end
    [z,i] = sort(x.(f)(is)); is = is(i);
end
ny = numel(y.(fm{1}));
js = 1:ny;
for k = 1:nm
    f = fm{k};
    if numel(y.(f))~=ny, error('Wrong field length.'), end
    [z,j] = sort(y.(f)(js)); js = js(j);
end

% define new fields in x

for k = 1:nc
    f = fc{k};
    if numel(x.(f))~=nx, error('Wrong field length.'), end
    x.([f 'm']) = NaN(size(x.(f)),class(x.(f)));
end

% find matching items

i = 1;
j = 1;
while i<=nx % for each item in x,
    ii = is(i);
    for jt = j:ny % look for a match in y
        jj = js(jt);
        for k = nm:-1:1
            f = fm{k};
            d = x.(f)(ii)-y.(f)(jj);
            if d~=0 % a field doesn't match
                break 
            end
        end
        if d==0 % all fields match
            for k = 1:nc
                f = fc{k};
                x.([f 'm'])(ii) = y.(f)(jj);
            end
            j = jt + 1;
            i = i + 1;
            break  % stop looking
        elseif d<0||jt==ny % no match exists for this item
            j = jt;
            i = i + 1;
            break  % stop looking
        end
   end
end
