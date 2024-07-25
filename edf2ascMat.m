function [status, msg] = edf2ascMat(fileName, exeFile, overwrite)
% edf2ascMat(fileName, exeFile, overwrite)
% convert edf file(s) to edf via system call
% fileName is either a string of one file name, or cell array of many
% exeFile is optional argument with full path to edf2asc.exe file
% overwrite [optional]: 1=overwrite the asc file, 0=never, []= ask each time
% default is "..\..\..\..\OneDrive - TCDUD.onmicrosoft.com\Other\edf_converter\@Edf2Mat\private\edf2asc.exe"

if ~exist('exeFile','var') || isempty(exeFile)
    exeFile = '"..\..\..\..\..\OneDrive - TCDUD.onmicrosoft.com\OtherOD\edf_converter\@Edf2Mat\private\edf2asc.exe"';
end
if ~exist('overwrite','var')
    overwrite = [];
end

if isa(fileName, 'char') % if one file
    
    [status, msg] = makeEdf(fileName, exeFile);
    
elseif isa(fileName, 'cell') % if a cell array of files
    
    for i = 1:length(fileName) % convert each cell separately
        [status, msg] = makeEdf(fileName{i}, exeFile, overwrite);
    end
    
else
    
    error('fileName must be a string or a cell array of strings')
    
end

end

function [status, msg] = makeEdf(fileName, exeFile, overwrite)
% makeEdf(fileName, exeFile)
% checks if edf file already exists, prompts for overwriting if so, creates
% if not
if ~exist(fileName,'file')
    error('edf file not found: %s', fileName)
end

ascName = [fileName(1:end-4) '.asc']; % name of asc file to check exists

if exist(ascName, 'file')
    
    if isempty(overwrite)
        % overwrite?
        r = input(sprintf('edf file %s already exists, overwrite? Y/N  ',ascName), 's');
        overwrite = strcmpi(r, 'y');
    end
    
    if overwrite
        
        fprintf('overwriting %s\n', ascName);
        system(['del ' ascName]); % delete old
        system([exeFile ' ' fileName]); % convert
        
    else
        fprintf('NOT overwriting %s\n', ascName);
    end
    
else % if no asc file already
    
    [status,msg] = system([exeFile ' ' fileName]);
    
end

end
        