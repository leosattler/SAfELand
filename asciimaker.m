function asciimaker(ModelFile, map, outputFolder, fileName)
%==========================================================================
% Function to generate ascii file of ginal susceptibility map.
%
% Inputs types: (string, array, string, string). 
% ModelFile = .tif file inside work folder, used to recuperate DEM infos. 
% These infos will be written on the header of the ascii output file as 
% follows:
% 1) NCOLS
% 2) NROWS
% 3) XLLCORNER
% 4) YLLCORNER
% 5) CELLSIZE
% 6) NODATA_VALUE (by default = -9999)
% Other inputs are:
% map = map of values (either from entire watershed or from scars, 
% consisting of q/t ratio values or shalstab classes)
% outputFolder = output folder location
% fileName = plot file name (suffix)
%==========================================================================
% Creating ascii output file name 
gridName = [outputFolder '\' fileName '.txt'];
%-------------------------------------------------------------------------- 
% Opening file
fid=fopen(gridName,'w');  
%-------------------------------------------------------------------------- 
% Retrieving headder info (from .tif model file)
info=geotiffinfo(ModelFile);  
%-------------------------------------------------------------------------- 
% Defining values for output file's header (from .tif model file)
[ny, nx] = size(map);  % NCOLS, NROWS
limits1=info.CornerCoords.X(end);  % XLLCORNER
limits2=info.CornerCoords.Y(end);  % YLLCORNER
cellsize=round(sqrt(info.FileSize/(info.Width*info.Height)));  % CELLSIZE
%-------------------------------------------------------------------------- 
% Writing header on output file
fprintf(fid, ...
    'ncols %d\nnrows %d\nxllcorner %16.4f\nyllcorner %16.4f\ncellsize %16.4f\nNODATA_value -9999\n', ... 
    nx, ny, limits1, limits2, cellsize);
%-------------------------------------------------------------------------- 
% Writing lines of data on output file
for i=1:ny
    fprintf(fid, ' %g',map(i,:));
    fprintf(fid, '\n');
end
%-------------------------------------------------------------------------- 
% Closing file 
fclose('all');
end