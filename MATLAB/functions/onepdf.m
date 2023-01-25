%%  Onepdf
%
% Create a single PDF document from multiple individual PDFs. Uses
% GPL-licensed PDF ToolKit (pdftk):
%       www.pdflabs.com
% 
% Requires pdftk.exe in any MATLAB path.
%       
% USAGE:
%   onepdf(pdfnames,outputname,deletefiles,Overwrite)
%
%   pdfnames: input filenames (can use wildcard *)
%   outputname: name of nultipage output PDF
%   deletefiles: removes original PDF files after combining
%   Overwrite: allow output file to be overwritten if it already exists

%---------------------------- 
% Dirk Albrecht 
% Version 1.0 14-Apr-2011 10:22:49
% Version 1.1 28-May-2021  allow filelist as dir structure (e.g.
%                           dir('*.pdf') for pdfnames argument
% Version 1.2 24-Jan-2023  changed D = dir(pdfnames) --> D = dir('*.pdf'). 
% 23:54:49                 Otherwise fails to get the PDFs full pathname, 
%                               error message: java.io.IOException: GetFullPathName failed 
%---------------------------- 

function onepdf(pdfnames,outputname,deletefiles,Overwrite)

if ~exist('pdfnames') || isempty(pdfnames) 
    [f,p] = uigetfile('MultiSelect','on');
    pdfnames = []; 
    for i=1:length(f); 
        pdfnames = cat(2,pdfnames,['"',p,char(f(i)),'" ']); 
    end  
end
if isstruct(pdfnames) && isfield(pdfnames,'name') % i.e. pdflist was generated from dir structure
    pdfnames = struct2cell(pdfnames);
    pdfnames = pdfnames(1,:);
    pdfnames = sprintf('"%s" ', pdfnames{:}); % now in the format of single string, '"filename1" "filename2" ...' 
end

if ~exist('outputname') || isempty(outputname) outputname = [p,'combinedPDF.pdf']; end
if ~exist('deletefiles') || isempty(deletefiles) deletefiles = true; end
if ~exist('Overwrite') || isempty(Overwrite) Overwrite = true; end

if isempty(strfind(outputname,'.pdf')) outputname = [outputname,'.pdf']; end

D = dir(outputname); OutputExist = size(D,1)>0;
if OutputExist
    if ~Overwrite
        error(['Output file ',PDFname,' exists.']);
        return;
    end
end

pdfpath = which('pdftk.exe');
if isempty(pdfpath)
    error('Can''t find pdftk.exe to create multipage pdf.  Place this file (http://www.accesspdf.com/pdftk/) in the matlab path.');
    return;
end

com = sprintf('!"%s" %s cat output %s',pdfpath,pdfnames,outputname);

currentpath = pwd;
if currentpath(1) == '\'
    p = path; lastpath = p(max(find(p == ';'))+1:length(p));
    cd(lastpath);
end

if (OutputExist & Overwrite) delete(outputname); end
eval(com);

disp(['Converted to PDF file: ',outputname]);
cd(currentpath);

D = dir(outputname); OutputExist = size(D,1)>0;
if OutputExist & deletefiles
    disp('Deleting single pdf pages');
    % Updated 2023-01-24 by VLK
        % D = dir(pdfnames); D = {D.name}';                  
    D = dir('*.pdf'); D = {D.name}';                    
    filestodelete = D(find(~strcmp(D,outputname)));     % exclude the combined pdf
    oldstate = recycle; recycle on;                     % turn on recycle bin, just in case
    for i=1:length(filestodelete)
        delete(char(filestodelete(i)));                 % delete pdfs
    end;
    recycle(oldstate);                                  % reset recycle bin state
end
end