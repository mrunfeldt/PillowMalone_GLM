%LOAD_RIPPLE loads entire contents of a DMR envelope file
%   RIPPLE = LOAD_RIPPLE(RIPPATH) loads the entire spectral-temporal
%   envelope of a .spr file and returns the results as the matrix RIPPLE
%   (rows = frequency samples, columns = time samples).
%
%   RIPPLE = LOAD_RIPPLE(RIPPATH,SPRTYPE) allows the user to specify the
%   data type in which the envelope data is stored (default: 'float').
%
%   RIPPLE = LOAD_RIPPLE(RIPPATH,SPRTYPE,TRIGGER,FS) allows the user to
%   provide a trigger vector that will enable the function to preallocate
%   memory for RIPPLE without checking the number of blocks in the file.
%
%   [RIPPLE,T_VEC] = LOAD_RIPPLE(...) returns a vector of timestamps (in
%   msec) for each column of the returned RIPPLE variable. 
%
%
%   Written by Jonathan Shih 1-2-2008
% Modified MJRunfeldt 2105_12_22 to output 'faxis' variable, which is
% loaded via the param path
function [ripple,t_vec,faxis] = load_ripple(rippath,sprtype,trigger,fs)

%Checking input parameters
if(~exist('sprtype','var'))
    sprtype = 'float';
end

%Loading parameter file
ind = findstr(rippath,'.spr');
parampath = [rippath(1:ind(1)-1) '_param.mat'];
load(parampath); % "faxis" is loaded from param file

%Finding length of file so that we can preallocate memory for RIPPLE
fid = fopen(rippath);
if(~exist('trigger','var'))
    blockCount = 0;
    while(~feof(fid))
        temp = fread(fid,NT*NF,sprtype);
        if(length(temp) == NT*NF)
            blockCount = blockCount + 1;
        end
    end
else
    blockCount = length(trigger);
    trig_space = 1000*mean(diff(trigger))/fs; %finding trigger spacing
end

%Initializing RIPPLE and T_VEC variable
ripple = zeros(NF,NT*blockCount);
if(exist('trigger','var'))
    t_vec = zeros(1,NT*blockCount);
else
    t_vec = [];
end

%Reading envelope file
frewind(fid);
fprintf('Loading ripple . . . ');
for block = 1:blockCount    
%     temp = fread(fid,NT*NF,sprtype);
%     if(length(temp) ~= length(ripple(((NT*NF*(block-1)+1)):(NT*NF*block))))
%         block
%         keyboard;
%     end
%     ripple(((NT*NF*(block-1)+1)):(NT*NF*block)) = temp;
    ripple(((NT*NF*(block-1)+1)):(NT*NF*block)) = fread(fid,NT*NF,sprtype);
    
    if(exist('trigger','var'))
        t_start = 1000*trigger(block)/fs;
        t_block = linspace(t_start,t_start+trig_space,NT+1);
        t_vec((NT*(block-1)+1):(NT*block)) = t_block(1:NT);
    end
    
%     %Updating progress
%     if(mod(block,200) == 0)
%         fprintf('Finished loading block %d of %d\n',block,blockCount);
%     end
end
fprintf('finished.\n');
fclose(fid);