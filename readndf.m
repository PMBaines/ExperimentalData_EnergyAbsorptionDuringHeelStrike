function [x,y,z] = readndf(name)

% Read 3D Optotrakdata.
    [fid,message]=fopen(name);

if fid == -1
    error(message);
end

% note: last parameter in following call is the letter 'l' for 'little
% endian' byte ordering

fid=fopen(name,'r','l');

% read header
filetype    = fread (fid, 1, 'char');
items       = fread (fid, 1, 'int16');
subitems    = fread (fid, 1, 'int16');
frames      = fread (fid, 1, 'int32');
freq        = fread (fid, 1, 'float32');
comment     = fread (fid, 60, 'char');
rest        = fread (fid, 183,'char');
is          = items*subitems;

% read data
data        = fread(fid, is*frames, 'float32');
index       = find(abs(data) > 3e+028);
data(index) = data(index).*NaN;
data        = reshape(data, is, frames)';

if subitems == 1
x = data(:, 1:1:is);
elseif subitems == 2
z = data(:, 1:2:is);
y = data(:, 2:2:is);
elseif subitems == 3
x = data(:, 1:3:is);
y = data(:, 2:3:is);
z = data(:, 3:3:is);
end

fclose(fid);