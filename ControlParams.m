%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       SSSSS    AAA    fff EEEEEEE LL                           dd 
%      SS       AAAAA  ff   EE      LL        aa aa nn nnn       dd 
%       SSSSS  AA   AA ffff EEEEE   LL       aa aaa nnn  nn  dddddd 
%           SS AAAAAAA ff   EE      LL      aa  aaa nn   nn dd   dd 
%       SSSSS  AA   AA ff   EEEEEEE LLLLLLL  aaa aa nn   nn  dddddd 
%
%          Soil Susceptibility Analysis for Estimating Landslides          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program status. At each parameter combination, the status of the run is
% saved. In the case of interruption before completion, user may return to
% last combination performed. Possible values are 0 (zero - program starts
% from begining) or 1 (one - program starts from interruption). If program
% was interrupted and RunInterrupted = 1, the inputs must be identical to
% the ones previously chosen.
%
RunInterrupted = 0; % program status (0 or 1)
%--------------------------------------------------------------------------
% Location for saving outputs. If output folder is not specified, the 
% default location is the same where program files are saved (pwd). If
% chosen folder already exists, it is EXCLUDED and recreated.
%
% Output folder:
%outputFolder=''; % specified
outputFolder=pwd; % same as program folder
%--------------------------------------------------------------------------
% Files of slope [rad] and upslope contributing area [m2] maps for 
% entire watershed (1) and for scars (2). The entire file path must be 
% specified. See mannual for instructions on how to correctly generate map 
% of scars from ArcGis, if needed.
%
% Entire basin:
thetamap1 = 'basin_slp.tif';
flowaccmap1 = 'basin_acc.tif ';
%
% Scars:
thetamap2 = 'scars_slp.tif';
flowaccmap2 = 'scars_acc.tif';
%
b=2.;
%--------------------------------------------------------------------------
% Geotechnical parameters for SHALSTAB calculation. They can be either
% lists of values or maps of spatially distributed values. Lists: numerical
% values separated by space; empty list is not accepted. Maps: .tif files;
% name of file must contain entire path. If both are (uncorrectly) defined
% before run, the program ignores list values and performs simulation using
% map only.
%
% phi = soil friction angles [degrees]
phiValues = [20 22 24 26 28 30 32 34 36 38 40 42 44 46];  %list of values
%phiValues = 'phi.tif';  %grid
%
% rhoS = soil density [kg/m3]
rhoSValues = [1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500];   %list of values                        
%rhoSValues = 'rhoS.tif';  %grid
%
% c = soil cohesion [Pa]
cValues = [0 1000 2000 3000 4000 5000 6000 7000 8000];  %list of values
%cValues = 'c.tif';  %grid
%
% z = soil thickness [m]
zValues = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];  %list of values
%zValues = 'z.tif';  %grid
%--------------------------------------------------------------------------
% Resolution = number of digits after decimal point to be truncated for 
% each pixel's q/T value, used to calculate number of possible classes and
% its intervals for validation function; the larger the resolution, greater
% is the number of defined classes. Default value is 6,
%
resolution = 6;
%--------------------------------------------------------------------------
% Number of outputs to be saved in Results folder from list of best 
% performing soil parameter combinations; if chosen number is greater than
% possible number of combinations, the top result is automatically chosen.
% Default value is 5.
%
NResults = 5;
%%
% Running the program: generating and saving outputs
outputgenerator(RunInterrupted, outputFolder, ...
                thetamap1, flowaccmap1, thetamap2, flowaccmap2, b,...
                phiValues, rhoSValues, cValues, zValues, ...
                resolution, NResults);
%%
% To reset all inputs
fclose('all');
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%