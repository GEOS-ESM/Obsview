function makedoc(varargin)

% MAKEDOC Update code-dependent parts of OBSTOOLS documentation
%
%   MAKEDOC
%   MAKEDOC(srcpath)
%   MAKEDOC(srcpath,htmldir) 
%
% Requires M2HTML (http://www.artefact.tk/software/matlab/m2html)
%
% Instructions:
% ------------
% - Install m2html and make sure it is in the Matlab search path 
%   (try: help m2html)
% - Check out the OBSTOOLS source code to a local directory
%   (for example, ./obsview')
% - Make sure this local directory is in the Matlab search path 
%   (try: which('obsview'))
% - Run this function
%   (try: makedoc, without any arguments)

curdir = pwd;  % current working directory
srcpath = fileparts(which('obsview')); % default source code directory
htmldir = 'obstools'; % default html output directory (relative to srcpath/..)

if nargin
   srcpath = varargin{1};
   if nargin>1
      htmldir = varargin{2};
   end
end

cd(srcpath)  % the code is here 
cd('..')     % this directory must be writeable

% create .m file documentation

[path,srcdir] = fileparts(srcpath);
m2html('mfiles',srcdir,'htmldir',[htmldir filesep 'code']);

% create kx, kt, qcx, qch tables

maketables(htmldir)

% go back where we started:

cd(curdir)

%-------------------------------------------------------------------

function maketables(dirname)

% This program creates:
%
%      kxtable.html
%      kttable.html
%      qchtable.html
%      qcxtable.html

[KXS,KTS,QCHS,QCXS] = dconfig('KXS','KTS','QCHS','QCXS');

% create kxtable.html
% -------------------

fname = [dirname filesep 'kxtable.html'];
disp(['Creating HTML file ' fname]) 
fid = fopen(fname,'w');

fprintf(fid,'%s','<html><head><title>Data Source Index Table</title>');
fprintf(fid,'%s','<link rel="stylesheet" type="text/css" href="./obstools.css">');
fprintf(fid,'%s','</head><body>');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><b>Data Source Index Table</b></caption>');
[w,i] = sort([KXS.value]);
for k = i,
   fprintf(fid,['<tr><td>%d</td><td>' KXS(k).id '</td>\n'],KXS(k).value);
end
fprintf(fid,'%s','</center></table>');


fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

% create kttable.html
% -------------------

fname = [dirname filesep 'kttable.html'];
disp(['Creating HTML file ' fname]) 
fid = fopen(fname,'w');

fprintf(fid,'%s','<html><head><title>Data Type Index Table</title>');
fprintf(fid,'%s','<link rel="stylesheet" type="text/css" href="./obstools.css">');
fprintf(fid,'%s','</head><body>');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><b>Data Type Index Table</b></caption>');
[w,i] = sort([KTS.value]);
for k = i,
   units = KTS(k).units;
   if strcmp(units,'%'), units = '%%'; end
   fprintf(fid,['<tr><td>%d</td><td>' KTS(k).id ' [' units ']</td>\n'],KTS(k).value);
end
fprintf(fid,'%s','</center></table>');

fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

% create qchtable.html
% -------------------

fname = [dirname filesep 'qchtable.html'];
disp(['Creating HTML file ' fname]) 
fid = fopen(fname,'w');

fprintf(fid,'%s','<html><head><title>Quality Control History Marks</title>');
fprintf(fid,'%s','<link rel="stylesheet" type="text/css" href="./obstools.css">');
fprintf(fid,'%s','</head><body>');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><b>Quality Control History Marks</b></caption>');

[w,i] = sort([QCHS.value]);
for k = i,
   fprintf(fid,['<tr><td>%d</td><td>' QCHS(k).id '</td>\n'],QCHS(k).value);
end
fprintf(fid,'%s','</center></table>');

fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);

% create qcxtable.html
% -------------------

fname = [dirname filesep 'qcxtable.html'];
disp(['Creating HTML file ' fname]) 
fid = fopen(fname,'w');

fprintf(fid,'%s','<html><head><title>Quality Control Exclusion Marks</title>');
fprintf(fid,'%s','<link rel="stylesheet" type="text/css" href="./obstools.css">');
fprintf(fid,'%s','</head><body>');

fprintf(fid,'%s','<center><table border>');
fprintf(fid,'%s','<caption><b>Quality Control Exclusion Marks</b></caption>');

[x,i] = sort([QCXS.value]);
for k = i,
   fprintf(fid,['<tr><td>%d</td><td>' QCXS(k).id '</td>\n'],QCXS(k).value);
end
fprintf(fid,'%s','</center></table>');

fprintf(fid,'%s\n','<p><hr>');
fprintf(fid,'%s\n',['<small><i>Generated by ' upper(mfilename) ' on ' datestr(now,0) '.']);
fprintf(fid,'%s\n','</body></html>');

fclose(fid);
