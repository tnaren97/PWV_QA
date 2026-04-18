%% Load Dat
% Author: Tarun Naren
% Date: October 3rd 2022
% Updated: April 4th 2024 by Tarun Naren
% Function to load in binary dat file data. 
function v = load_dat(filename, res)
    arguments
        filename (1,1) string
        res (1,:) double          % matrix size (specify as 1 by n vector, Ex: [256 256 256])
    end

    try
        m = memmapfile(filename, 'Format', {'int16', prod(res), 'data'});
        v = single(reshape(m.Data.data, res));
        v = flip(flip(v,1),3); % orient in LAS+ direction
    catch ME
        warning('Error loading file %s: %s', filename, ME.message);
    end
end 