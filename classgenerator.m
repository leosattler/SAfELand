function map = classgenerator(THETA, FLOWACC, phi, rhoS, c, z)
%==========================================================================
% Function to generate class susceptibility map.
%
% Input types: (array, array, array, array, array, array). 
% THETA = slope map [rad]
% FLOWACC = contributing are map [m2]
% phi = friction angle of soil [graus]
% rhoS = soil density [kg/m3]
% c = soil cohesion [Pa]
% z = soil thickness [m]
%==========================================================================
% Defining variables for entire DEM
rhoagua=1000;    % water density [kg/m3]
g=9.8;           % gravity [m/s2]
%--------------------------------------------------------------------------
% Creating grid of susceptibility map
map=zeros(size(THETA));  % grid for all classes
%--------------------------------------------------------------------------
% Initiating loop (only pixels inside DEM, avoiding NoData)
dem=find(THETA >= 0);
for i=1:length(dem)
    pixel=dem(i);
    % Constants C1 and C2 for cohesion fraction 
    C1=c(pixel)/((cos(THETA(pixel))^2)*rhoS(pixel)*z(pixel)*g); 
    C2=c(pixel)/((cos(THETA(pixel))^2)*tan(phi(pixel))*z(pixel)*rhoagua*g); 
    % Finding pixels under ALWAYS STABLE condition   
    if tan(THETA(pixel)) >= (tan(phi(pixel)) + C1)
        map(pixel)=-10; 
    % Finding pixels under ALWAYS UNSTABLE condition    
    elseif tan(THETA(pixel)) <= (tan(phi(pixel)) * (1-(rhoagua/rhoS(pixel))) + C1)
        map(pixel)=10;
    % Avoiding values=log(0) and finding INTERMEDIARY conditions    
        elseif (1-(tan(THETA(pixel))/tan(phi(pixel)))) > 0 || C2 > 0
        % Temporary values to find intermediary class values
        tmp=(1-(tan(THETA(pixel))/tan(phi(pixel))));
        % Classifying
        map(pixel)=log10( sin(THETA(pixel)) * ((rhoS(pixel)/rhoagua) * ... 
        tmp + C2) / FLOWACC(pixel) );
    end
end
%--------------------------------------------------------------------------
% Attributing value -9999 to pixels outside watershed (NoData code)
map(THETA<0)=-9999;   
end