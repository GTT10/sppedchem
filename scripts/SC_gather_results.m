%% This subroutine gathers results from a SpeedCHEM's output file, 
%% output.dat, or output.bin in case a binary version is present
%% Federico Perini, 2013

% INPUT: folder = full or relative folder path
%        (optional) = filename in case it is not output

function [solution] = SC_gather_results(folder,varargin)

% Choose desired filename in case given as input
if (length(varargin)>0),
    filename = varargin{1};
else
    filename = 'output';
end

% Initialize output structure
solution.t = [];
solution.T = [];
solution.Y = [];

% Look for output.dat file in the folder
use_binary    = 1;
file_pathname = [folder '\' filename '.bin'];
file_path     = dir(file_pathname);

if (length(file_path)<1),
    use_binary    = 0;
    file_pathname = [folder '\' filename '.dat'];
    file_path     = dir(file_pathname);    
    if (length(file_path)<1),
    error([' No SpeedCHEM output found in folder ' folder]);
    end
end

% Read SpeedCHEM output in binary format
if (use_binary),

    % First of all, open file
    fid = fopen(file_pathname,'r');
    
    % Second, read output file dimensions: case_no, points_per_case,
    % nspecies
    [dummy, count]  = fread(fid,4);
    values          = fread(fid,3,'int');
    case_no         = values(1);
    points_per_case = values(2);
    nspecies        = values(3);
    [dummy, count]  = fread(fid,4);
    
    % Read all case information
    [dummy, count]  = fread(fid,4);
    values          = fread(fid,case_no*points_per_case*(nspecies+2),'double');
    [dummy, count]  = fread(fid,4);
    
    % Reshape values
    values = reshape(values, nspecies+2, points_per_case, case_no);
    
    for i = 1:case_no,
        solution(i).t = squeeze(values(1,:,i));
        solution(i).T = squeeze(values(2,:,i));
        solution(i).Y = squeeze(values(3:nspecies+2,:,i))';
    end
    
    % Close file
    fclose(fid);
    
    return
    
else
    
% Open output.dat for read. Structured read not possible.
fid = fopen(file_pathname,'r');

% Read while possible
nline = 0;
while ~feof(fid),
    
    line    = fgets(fid);
    
    case_no = sscanf(line,' Case  %d');

    if (~isempty(case_no)),
        nline = 0;
        current_case = case_no;
        
    else
        
        values = sscanf(line,'%f');
        nline  = nline + 1;
        
        solution(current_case).t(nline) = values(1);
        solution(current_case).T(nline) = values(2);
        solution(current_case).Y(nline,:) = values(3:end);
        
    end
    
end


% Close file output.dat
fclose(fid);

end

return