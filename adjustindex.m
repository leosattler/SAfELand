function ai = adjustindex(mapBasin, mapScars)
%==========================================================================
% Auxiliar function for the suceptibility map validation based on the
% adjustment index (AI) analysis: 
% AI = %PHS/PS, where:
% PHS (Pixels with Highest Suceptibility) = x% of watershed pixels that 
% both have the highest suceptibility values and reside inside a scar. 
% By default, this validation algorithm considers 4 values of x:
% x=5%, 10%, 20% and 30%.
% PS (Pixels inside Scars) = Amount of all pixels belonging to scars.
% This routine returns a list of all adjustment indexe values with the
% considered x values.
%
% Input types: (array, array).
% mapBasin = Basin's suceptibility map 
% mapScars = Scars' suceptibility map 
%==========================================================================
% Getting DEM pixel values (avoiding NoData)
mapBasin=mapBasin(mapBasin~=-9999);
%-------------------------------------------------------------------------- 
% Sorting pixels from least susceptible to most susceptible (indexes only)
[~,indexes]=sort(mapBasin);
%-------------------------------------------------------------------------- 
% Generating list of percenteages of pixels inside watershed
percentages=[round((5/100)*length(mapBasin)), round((10/100)*length(mapBasin)), ...
    round((20/100)*length(mapBasin)), round((30/100)*length(mapBasin))];
%-------------------------------------------------------------------------- 
% Generating (yet empty) list of adjustment indexes (AIs)
ai=zeros(1,length(percentages));
%-------------------------------------------------------------------------- 
% Getting Scar pixel values (avoiding NoData; indexes only)
PC=find(mapScars~=-9999);
%-------------------------------------------------------------------------- 
% Looping over x% values and calculating the resulting AI for each
for i=1:length(percentages)
    % Defining x% percentage for the loop
    pctg = percentages(i);
    % Finding intersection between x% of most basin's susceptible pixels 
    % and pixels belonging to scars
    aiLoop=length(intersect(indexes(1:pctg),PC));
    % Saving AI value inside AI list
    ai(i) = 100*(aiLoop/length(PC));
end
end