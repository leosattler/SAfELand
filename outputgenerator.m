function outputgenerator(RunInterrupted, outputFolder, thetamap1, flowaccmap1, thetamap2,... 
    flowaccmap2, b, phiArg, rhoSArg, cArg, zArg, resolution, NResults)
%==========================================================================
% Central function. Here the program receives inputs and calls each
% function.
%
% Input types: (double, string, string, string, string, string, 
%   string ou array, string ou array, string ou array, string ou array,
%   double, double)
% RunInterrupted = program status (run or re-run)
% outputFolder = complete path of putput folder
% thetamap1 = file name for watershed's slope [rad]
% flowaccmap1 = file name for watershed's upslope contributing area [m2]
% thetamap2 = file name for scars' slope [rad]
% flowaccmap2 = file name for scars' upslope contributing area [m2]
% phiArg = soil friction angles [degrees]
% rhoSArg = soil density [kg/m3]
% cArg = soil cohesion [Pa]
% zArg = soil thickness [m]
% resolution = resolution of class interval for validation function
% NResults = number of combination results to be outputed 
%==========================================================================
% Turning off matlab unecessary messages
% how to retrieve messages' code: [msg,msgID] = lastwarn; 
%--------------------------------------------------------------------------
warning('OFF', 'map:geotiff:undefinedGTModelTypeGeoKey'); % geotiff model
warning('OFF', 'MATLAB:table:ModifiedVarnames'); % table header 
warning('OFF', 'MATLAB:xlswrite:AddSheet'); % xls output
warning('OFF', 'MATLAB:singularMatrix'); % matrix precision
warning('OFF', 'MATLAB:legend:IgnoringExtraEntries'); % extra legend entries
warning('OFF', 'MATLAB:MKDIR:DirectoryExists'); % existing directory
%==========================================================================
% Checking if inputs RunInterrupted, outputflder, resolution e NResults are
% adequate for the progeam (if are positive numbers and if specified
% folders exist)
%--------------------------------------------------------------------------
% Returning if RunInterrupted is not adequate: 
% (negative number or non-number or non-integer)
if RunInterrupted<0 || ~isnumeric(RunInterrupted) || (floor(RunInterrupted) ~= RunInterrupted)
    clc;
    fprintf('Warning: RunInterrupted provided is invalid.')
    fprintf('\n');
    fprintf('Only non-negative integer values are accepted.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return
end
%--------------------------------------------------------------------------
% Returning if output folder is non-existent
if ~isdir(outputFolder)  
    clc;
    fprintf('Warning: the output folder provided does not exist.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return
end
%--------------------------------------------------------------------------
% Returning if resolution is not adequate:
% (negative number or non-number or non-integer)
if resolution<0 || ~isnumeric(resolution) || (floor(resolution) ~= resolution)
    clc;
    fprintf('Warning: resolution provided is invalid.')
    fprintf('\n');
    fprintf('Only non-negative integer values are accepted.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return
end
%--------------------------------------------------------------------------
% Returning if NResults is not adequate:
% (negative number or non-number or non-integer)
if NResults<0 || ~isnumeric(NResults) || (floor(NResults) ~= NResults)
    clc;
    fprintf('Warning: NResults provided is invalid.')
    fprintf('\n');
    fprintf('Only non-negative integer values are accepted.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return
end
%==========================================================================
% Begining message of program
fprintf('Importing files...');
fprintf('\n');
pause(1.);
%==========================================================================
% Checking if files exist and importing grids if so
%--------------------------------------------------------------------------
% Watershed slope
if exist(thetamap1,'file') == 0 % Checando se arquivo existe
    clc;
    fprintf('Warning: the file for thetamap1 was not found.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return;
else
    THETA1 = imread(thetamap1);
end
%--------------------------------------------------------------------------
% Watershed upslope contributing area
if exist(flowaccmap1,'file') == 0 % Checando se arquivo existe
    clc;
    fprintf('Warning: the file for flowaccmap1 was not found.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return;
else
    FLOWACC1 = imread(flowaccmap1);
    FLOWACC1 = FLOWACC1/b;
end
%--------------------------------------------------------------------------
% Scars slope
if exist(thetamap2,'file') == 0 % Checando se arquivo existe
    clc;
    fprintf('Warning: the file for thetamap2 was not found.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return;
else
    THETA2 = imread(thetamap2);
end
%--------------------------------------------------------------------------
% Scars upslope contributing area
if exist(flowaccmap2,'file') == 0 % Checando se arquivo existe
    clc;
    fprintf('Warning: the file for flowaccmap2 was not found.')
    fprintf('\n');
    fprintf('Exited without calculation.')
    fprintf('\n');    
    return;
else
    FLOWACC2 = imread(flowaccmap2);
    FLOWACC2 = FLOWACC2/b;
end
%==========================================================================
% Checking if inputs are string (name of .tif file) or lists
%--------------------------------------------------------------------------
% phi
if ischar(phiArg)
    % Checking if file for phi grid exists
    if exist(phiArg,'file') == 0
        clc;
        fprintf('Warning: the file for phi was not found.')
        fprintf('\n');
        fprintf('Exited without calculation.')
        fprintf('\n');    
        return;
    else
        phiValues = imread(phiArg); % if .tif file
    end        
else
    phiValues = phiArg; % if list
end
%--------------------------------------------------------------------------
% rhoS
if ischar(rhoSArg)
    % Checking if file for rhoS grid exists
    if exist(rhoSArg,'file') == 0
        clc;
        fprintf('Warning: the file for rhoS was not found.')
        fprintf('\n');
        fprintf('Exited without calculation.')
        fprintf('\n');   
        return;
    else
        rhoSValues = imread(rhoSArg); % if .tif file
    end
else
    rhoSValues = rhoSArg; % if list
end
%--------------------------------------------------------------------------
% c
if ischar(cArg)
    % Checking if file for c grid exists
    if exist(cArg,'file') == 0
        clc;
        fprintf('Warning: the file for c was not found.')
        fprintf('\n');
        fprintf('Exited without calculation.')
        fprintf('\n');    
        return;
    else 
        cValues = imread(cArg); % if .tif file
    end
else
    cValues = cArg; % if list
end
%--------------------------------------------------------------------------
% z
if ischar(zArg)
    % Checking if file for z grid exists
    if exist(zArg,'file') == 0
        clc;
        fprintf('Warning: the file for z was not found.')
        fprintf('\n');
        fprintf('Exited without calculation.')
        fprintf('\n');    
        return;
    else
        zValues = imread(zArg); % if .tif file
        % Correcting possible values of z=0 and returnin 0.1,
        % since z is the denominator in the equation for q/T.      
        zValues(zValues==0) = 0.1;
    end
else
    zValues = zArg; % if list
end
%==========================================================================
% Working the inputs. If n# of rows=1, input is list. If >1, input is grid 
% (and won't be selected for cominations).
%--------------------------------------------------------------------------
[rowsphi, colsphi] = size(phiValues);
[rowsrhoS, colsrhoS] = size(rhoSValues);
[rowsc, colsc] = size(cValues);
[rowsz, colsz] = size(zValues);
%==========================================================================
% Defining limits of for loop of each input. If input is grid, its end-of-loop 
% value will be 1 (is not part of loop so won't be used for combination). 
% If input is list, its end-of-loop value is the size of list. 
% Also, auxiliar variables are generated (strcheck) for checking log.bin file. This
% file contains informations of each input parameter (if list or gid, and
% the values in case of list).
%--------------------------------------------------------------------------
% phi
if rowsphi~=1
    ifinal=1; % if grid
    phibin = 0;
else
    ifinal=colsphi; % if list
    phibin = colsphi;
end
%--------------------------------------------------------------------------
% rhoS
if rowsrhoS~=1
    jfinal=1; % if grid
    rhoSbin=0;
else
    jfinal=colsrhoS; % if list
    rhoSbin=colsrhoS;
end
%--------------------------------------------------------------------------
% c
if rowsc~=1
    kfinal=1; % if grid
    cbin=0;
else
    kfinal=colsc; % if list
    cbin=colsc;
end
%--------------------------------------------------------------------------
% z
if rowsz~=1
    lfinal=1; % if grid
    zbin=0;
else
    lfinal=colsz; % if list
    zbin=colsz;
end
%==========================================================================
% Organizing folders and checking variables of program start up: deleting 
% and recreating output folders , observing consistency of inputs, checking
% log and etc.
%--------------------------------------------------------------------------
% Naming folders for output: integral and adjustiment index methods +
% Result
outputFolderArea = [outputFolder '\Percentage_Area'];
outputFolderAdjIndx = [outputFolder '\Adjustment_Index'];
outputFolderResults = [outputFolder '\Results'];
%==========================================================================
% Checking run stage. If it starts from zero, existence of output folders
% is checked and deleted if one of same name exists. If run was
% interrupted, existence of binary files is checked along with the
% consistency between last run inputs and actual run inputs. Output folders
% are created (if do not exist).
%--------------------------------------------------------------------------
% If run starts from zero
if ~RunInterrupted 
    %----------------------------------------------------------------------
    % Deleting only folders created by the program (if they exist)    
    if isdir(outputFolderArea) 
        rmdir(outputFolderArea, 's'); 
    end 
    %
    if isdir(outputFolderAdjIndx) 
        rmdir(outputFolderAdjIndx, 's'); 
    end
    %
    if isdir(outputFolderResults) 
        rmdir(outputFolderResults, 's'); 
    end
    %----------------------------------------------------------------------
    % Deleting only files created by the program (if they exist)   
    if exist([outputFolder '\log.txt'], 'file') 
        delete([outputFolder '\log.txt']);
    end
    if exist([outputFolder '\log.bin'], 'file')
        delete([outputFolder '\log.bin']);
    end
    if exist([outputFolder '\Percent_Area_outputs.bin'], 'file')
        delete([outputFolder '\Percent_Area_outputs.bin']);
    end
    if exist([outputFolder '\Adjust_Index_outputs.bin'], 'file')
        delete([outputFolder '\Adjust_Index_outputs.bin']);
    end
    %----------------------------------------------------------------------
    % Creating output folders
    % Folders of each validation method
    mkdir(outputFolderArea);
    mkdir(outputFolderAdjIndx);   
    mkdir(outputFolderResults);
    % Integral method: Minimum and mode 
    mkdir([outputFolderArea '\Scars_Min']);  % minimo 
    mkdir([outputFolderArea '\Scars_Mode']); % moda
    %mkdir([outputFolderArea '\Scars_Raw']);  % raw 
    % To save ascii grids 
    mkdir([outputFolderArea '\Scars_Min\Grids']); 
    mkdir([outputFolderArea '\Scars_Mode\Grids']);
    %mkdir([outputFolderArea '\Scars_Raw\Grids']); 
    % To save figures of susceptbility maps 
    mkdir([outputFolderArea '\Scars_Min\Sucept_Plots']); 
    mkdir([outputFolderArea '\Scars_Mode\Sucept_Plots']);
    %mkdir([outputFolderArea '\Scars_Raw\Sucept_Plots']); 
    % To save validation plots
    mkdir([outputFolderArea '\Scars_Min\Validation']);  
    mkdir([outputFolderArea '\Scars_Mode\Validation']); 
    %mkdir([outputFolderArea '\Scars_Raw\Validation']);  
    % To save ascii grids 
    mkdir([outputFolderAdjIndx '\Grids']); 
    % To save figures of susceptibility maps 
    mkdir([outputFolderAdjIndx '\Sucept_Plots']);
    % Para guardar plots de validacao
    mkdir([outputFolderAdjIndx '\Validation']); 
    % Rankings of 5%, 10%, 20% e 30% for Adjustment Index method
    mkdir([outputFolderAdjIndx '\Validation\Ranking_5%']); % 5% 
    mkdir([outputFolderAdjIndx '\Validation\Ranking_10%']); % 10%
    mkdir([outputFolderAdjIndx '\Validation\Ranking_20%']); % 20%
    mkdir([outputFolderAdjIndx '\Validation\Ranking_30%']); % 30%
    %----------------------------------------------------------------------
    % Creating log.bin file to save inputs and run state (i.e. counter
    % value)
    fclose('all');    
    fileBinLog = fopen([outputFolder '\log.bin'],'w'); % Opening file
    if ~ischar(phiArg) 
        fwrite(fileBinLog, phiArg,'double');   % phi
    end
    if ~ischar(rhoSArg)
        fwrite(fileBinLog, rhoSArg,'double');  % rhoS
    end
    if ~ischar(cArg)
        fwrite(fileBinLog, cArg,'double');     % c
    end
    if ~ischar(zArg)
        fwrite(fileBinLog, zArg,'double');     % z
    end
    fwrite(fileBinLog, resolution, 'double');  % Resolution
    fwrite(fileBinLog, 0, 'double');           % Counter    
    fclose('all');  
    % Creating vector for update of counter in log.bin file during loop
    fileBinLog = fopen([outputFolder '\log.bin'],'r');
	fileBinColumn = fread(fileBinLog,'double');
    fclose('all');  
    %----------------------------------------------------------------------
    % Creating log.txt via information table for visualization of inputs
    % and run's current stage (same infos as log.bin)
    logTable = table(...
    {'Run Reinitiated at';'-';'phi';'rhoS';'c';'z';'-'; ... 
    'thetamap1';'flowaccmap1';'thetamap2';'flowaccmap2';'-'; ...
    'Resolution';'Counter';'-';'Run finished at'}, ...
    {datestr(datetime('now'));'-';mat2str(phiArg);mat2str(rhoSArg); mat2str(cArg);mat2str(zArg);'-'; ...
    thetamap1;flowaccmap1;thetamap2;flowaccmap2;'-'; ...
    num2str(resolution);'0';'-';'-'});
    logTable.Properties.VariableNames={'Log' 'File'};
    writetable(logTable,[outputFolder '\log.txt'],'Delimiter',' ');
    fclose('all'); 
%--------------------------------------------------------------------------
% if run was interrupted
else
    %----------------------------------------------------------------------
    % If folder is empty and without log.bin file (user incorrectly assigned      
    % RunInterrupted = 1 or chose wrong folder)
    if ~exist([outputFolder '\log.bin'], 'file')% Checking if log.bin file exists
        clc;
        fprintf('Warning: the output folder provided does not contain a log file.')
        fprintf('\n');
        fprintf('Check the path or reset RunInterrupted to 0.')
        fprintf('\n');   
        fprintf('Exited without calculation.')
        fprintf('\n');    
        return        
    %----------------------------------------------------------------------
    % Checking (via log.bin file) if inputs of current run are same as
    % last. If yes, run proceeds. If no, returns error message.
    %----------------------------------------------------------------------
    % Procedure to check if current list is same as last: each element of 
    % fileBinColumn (a column vector from log.bin file comprised of all
    % inputs of last run) with each list of inputs from current run. The
    % variable 'name of parameter+bin' (e.g., phibin) counts how many
    % entries there are in each parameters list (e.g., phiValues=[20 25 30] 
    % yields phibin=3). This procedure also checks if particular parameter 
    % was list and is now grid (and vice versa) using the log.bin file 
    % (via strcheck variable).
    else
        fclose('all');
        % Importing file (single column vector)
        fileBinLog = fopen([outputFolder '\log.bin'],'r');
        % Reading it. This vector will be used to update counter along loop
        fileBinColumn = fread(fileBinLog,'double');
        fclose('all');
        % Creating auxiliar variable to check consistency among runs
        if fileBinColumn(end) == -9999
            strcheck = length(fileBinColumn)-3;
        else
            strcheck = length(fileBinColumn)-2;
        end
        %------------------------------------------------------------------
        % Checking phi input
        k=1; 
        if ~ischar(phiArg) % If list, compare each entry
            for i=1:phibin
                if fileBinColumn(i) ~= phiArg(k) % If entry is different
                    % Run is cancled with error message
                    clc;
                    fprintf('Warning: current input for phi is not consistent with previous run.');
                    fprintf('\n');
                    fprintf('Please restart run from zero or change your inputs.');
                    fprintf('\n');
                    fprintf('Exited without calculation.')
                    fprintf('\n');   
                    % End of message and return
                    return;
                end
                k=k+1;
            end
        else % If grid, check if last run was also grid
            if (rhoSbin+cbin+zbin)==strcheck 
                % If yes, print attention (for file path) message 
                fprintf('Make sure input file for phi is same as before!');
                fprintf('\n');                
                pause(2);
                fprintf('Continuing... ');
                fprintf('\n');
                pause(1);
            else
                % If no, run is canceled with error message
                clc;
                fprintf('Warning: current input for phi is not consistent with previous run.');
                fprintf('\n');
                fprintf('Please restart run from zero or change your inputs.');
                fprintf('\n');
                fprintf('Exited without calculation.')
                fprintf('\n');   
                % End of message and return
                return;
            end
        end
        %------------------------------------------------------------------
        % Checking shoS input
        k=1;
        if ~ischar(rhoSArg) % If list, compare each entry
            for i=(1+phibin):(phibin+rhoSbin)
                if fileBinColumn(i) ~= rhoSArg(k) % If entry is different
                    % Run is cancled with error message
                    clc;
                    fprintf('Warning: current input for rhoS is not consistent with previous run.');
                    fprintf('\n');
                    fprintf('Please restart run from zero or change your inputs.');
                    fprintf('\n');
                    fprintf('Exited without calculation.')
                    fprintf('\n');   
                    % Fim da mensagem e retorno
                    return;
                end
                k=k+1;
            end
        else % If grid, check if last run was also grid
            if (phibin+cbin+zbin)==strcheck
                % If yes, print attention (for file path) message           
                fprintf('Make sure input file for rhoS is same as before!');
                fprintf('\n');                
                pause(2);
                fprintf('Continuing... ');
                fprintf('\n');
                pause(1);
            else
                % If no, run is canceled with error message
                clc;
                fprintf('Warning: current input for rhoS is not consistent with previous run.');
                fprintf('\n');
                fprintf('Please restart run from zero or change your inputs.');
                fprintf('\n');
                fprintf('Exited without calculation.')
                fprintf('\n');   
                % End of message and return
                return;
            end
        end
        %------------------------------------------------------------------
        % Checking c input
        k=1;
        if ~ischar(cArg) % % If list, compare each entry
            for i=(1+phibin+rhoSbin):(phibin+rhoSbin+cbin)
                if fileBinColumn(i) ~= cArg(k) % If entry is different
                    % Run is cancled with error message
                    clc;
                    fprintf('Warning: current input for c is not consistent with previous run.');
                    fprintf('\n');
                    fprintf('Please restart run from zero or change your inputs.');
                    fprintf('\n');
                    fprintf('Exited without calculation.')
                    fprintf('\n');   
                    % End of message and return
                    return;
                end
                k=k+1;
            end
        else % If grid, check if last run was also grid
            if (phibin+rhoSbin+zbin)==strcheck
                % If yes, print attention (for file path) message         
                fprintf('Make sure input file for c is same as before!');
                fprintf('\n');                
                pause(2);
                fprintf('Continuing... ');
                fprintf('\n');
                pause(1);
            else
                % If no, run is canceled with error message
                clc;
                fprintf('Warning: current input for c is not consistent with previous run.');
                fprintf('\n');
                fprintf('Please restart run from zero or change your inputs.');
                fprintf('\n');
                fprintf('Exited without calculation.')
                fprintf('\n');   
                % End of message and return
                return;
            end
        end
        %------------------------------------------------------------------
        % Checking z input
        k=1;
        if ~ischar(zArg) % % If list, compare each entry
            for i=(1+phibin+rhoSbin+cbin):(phibin+rhoSbin+cbin+zbin)
                if fileBinColumn(i) ~= zArg(k) % If entry is different
                    % Run is cancled with error message
                    clc;
                    fprintf('Warning: current input for z is not consistent with previous run.');
                    fprintf('\n');
                    fprintf('Please restart run from zero or change your inputs.');
                    fprintf('\n');
                    fprintf('Exited without calculation.')
                    fprintf('\n');   
                    % End of message and return
                    return;
                end
                k=k+1;
            end
        else % If grid, check if last run was also grid
            if (phibin+rhoSbin+cbin)==strcheck
                % If yes, print attention (for file path) message         
                fprintf('Make sure input file for z is same as before!');
                fprintf('\n');                
                pause(2);
                fprintf('Continuing... ');
                fprintf('\n');
                pause(1);
            else
                % If no, run is canceled with error message
                clc;
                fprintf('Warning: current input for z is not consistent with previous run.');
                fprintf('\n');
                fprintf('Please restart run from zero or change your inputs.');
                fprintf('\n');
                fprintf('Exited without calculation.')
                fprintf('\n');   
                % End of message and return
                return;
            end
        end
        %------------------------------------------------------------------
        % Checking if chosen resolution is same as last run. This is only
        % relevant if run is not complete (meaning last line at log.bin file
        % is not -9999). If run is complete, warning message is returned
        % and run proceeds with original resolution (as present at log.bin
        % file).
        if resolution ~= fileBinColumn(end-1) && fileBinColumn(end) ~= -9999
            clc;
            fprintf('Warning: current resolution is not consistent with previous run.');
            fprintf('\n');
            fprintf('Please restart run from zero or change your inputs.');            
            fprintf('\n');
            fprintf('Exited without calculation.')
            fprintf('\n');   
            return;
        elseif resolution ~= fileBinColumn(end-2) && fileBinColumn(end) == -9999
            clc;
            fprintf('The input resolution is not consistent with previous run.');
            fprintf('\n');
            fprintf(['Original resolution was ' num2str(fileBinColumn(end-2)) ' and will be used instead.']);
            fprintf('\n');
            pause(4.);
            resolution = fileBinColumn(end-2);
        end
        %------------------------------------------------------------------
        % Retrieving last counter (run stage before interruption)    
        if fileBinColumn(end) == -9999 % If run was completed
            LastCounter = fileBinColumn(end-1);
        else                          % If run was not completed
            LastCounter = fileBinColumn(end);
        end
        %------------------------------------------------------------------
        % Checking if log.txt file exists. If so, it is deleted and 
        % recreated. File will be rewritten from scratch with new LasCounter.   
        if exist([outputFolder '\log.txt'], 'file') % If exists
            % Deleting
            delete([outputFolder '\log.txt']); 
        end
        % If counter is ongoing (run not complete)
        if fileBinColumn(end) ~= -9999
            % Creating table
            logTable = table(...
            {'Run Reinitiated at';'-';'phi';'rhoS';'c';'z';'-'; ... 
            'thetamap1';'flowaccmap1';'thetamap2';'flowaccmap2';'-'; ...
            'Resolution';'Counter';'-';'Run finished at'}, ...
            {datestr(datetime('now'));'-';mat2str(phiArg);mat2str(rhoSArg); mat2str(cArg);mat2str(zArg);'-'; ...
            thetamap1;flowaccmap1;thetamap2;flowaccmap2;'-'; ...
            num2str(resolution);num2str(LastCounter);'-';'-'});
            logTable.Properties.VariableNames={'Log' 'File'};
            % Writing
            writetable(logTable,[outputFolder '\log.txt'],'Delimiter',' ');
        % If run is marked as interrupted and last line at log.bin file is -9999.
        % E.g. if user has only binary files and wishes to reproduce
        % outputs from last run (already completed).
        else 
            % Creating table
            logTable = table(...
            {'Run Reinitiated at';'-';'phi';'rhoS';'c';'z';'-'; ... 
            'thetamap1';'flowaccmap1';'thetamap2';'flowaccmap2';'-'; ...
            'Resolution';'Counter';'-';'Run finished at'}, ...
            {datestr(datetime('now'));'-';mat2str(phiArg);mat2str(rhoSArg); mat2str(cArg);mat2str(zArg);'-'; ...
            thetamap1;flowaccmap1;thetamap2;flowaccmap2;'-'; ...
            num2str(resolution);'complete';'-';'-'});
            logTable.Properties.VariableNames={'Log' 'File'};
            % Writing
            writetable(logTable,[outputFolder '\log.txt'],'Delimiter',' ');            
        end
        %------------------------------------------------------------------
        % Recreating folders if they were excluded (nothing is deleted!)       
        % Percentage Area folder
        if ~isdir(outputFolderArea)             
            mkdir(outputFolderArea);
        end
        % Scar min analysis
        if ~isdir([outputFolderArea '\Scars_Min'])
            mkdir([outputFolderArea '\Scars_Min']);  
        end
        % Scar mode analysis
        if ~isdir([outputFolderArea '\Scars_Mode'])
            mkdir([outputFolderArea '\Scars_Mode']); 
        end
        % Scar raw analysis (not implemented)
        %if ~isdir([outputFolderArea '\Scars_Raw'])
        %    mkdir([outputFolderArea '\Scars_Raw']);  
        %end         
        % To save ascii grids (Min)
        if ~isdir([outputFolderArea '\Scars_Min\Grids'])
            mkdir([outputFolderArea '\Scars_Min\Grids']); 
        end
        % To save ascii grids (Mode)
        if ~isdir([outputFolderArea '\Scars_Mode\Grids'])
            mkdir([outputFolderArea '\Scars_Mode\Grids']);
        end
        % To save ascii grids (Raw) (not implemented)
        %if ~isdir([outputFolderArea '\Scars_Raw\Grids'])
        %    mkdir([outputFolderArea '\Scars_Raw\Grids']);
        %end         
        % To save suscept. map figures via Shalstab model (Min)
        if ~isdir([outputFolderArea '\Scars_Min\Sucept_Plots'])
            mkdir([outputFolderArea '\Scars_Min\Sucept_Plots']); 
        end
        % To save suscept. map figures via Shalstab model (Mode)
        if ~isdir([outputFolderArea '\Scars_Mode\Sucept_Plots'])
            mkdir([outputFolderArea '\Scars_Mode\Sucept_Plots']);
        end
        % To save suscept. map figures via Shalstab model  (Raw) (not implemented)
        %if ~isdir([outputFolderArea '\Scars_Raw\Sucept_Plots'])
        %    mkdir([outputFolderArea '\Scars_Raw\Sucept_Plots']);
        %end
        % To save validation plots (Min)
        if ~isdir([outputFolderArea '\Scars_Min\Validation'])
            mkdir([outputFolderArea '\Scars_Min\Validation']);  
        end
        % To save validation plots (Mode)
        if ~isdir([outputFolderArea '\Scars_Mode\Validation'])
            mkdir([outputFolderArea '\Scars_Mode\Validation']); 
        end
        % To save validation plots (Raw) (not implemented)
        %if ~isdir([outputFolderArea '\Scars_Raw\Validation'])
        %    mkdir([outputFolderArea '\Scars_Raw\Validation']);                       
        %end 
        % Adjustment Index folder
        if ~isdir(outputFolderAdjIndx)                                   
            mkdir(outputFolderAdjIndx);
        end
        % To save ascii grids
        if ~isdir(([outputFolderAdjIndx '\Grids']))
            mkdir([outputFolderAdjIndx '\Grids']); 
        end
        % To save suscept. map figures via Shalstab model
        if ~isdir([outputFolderAdjIndx '\Sucept_Plots'])
            mkdir([outputFolderAdjIndx '\Sucept_Plots']);
        end
        % To save validation plots
        if ~isdir([outputFolderAdjIndx '\Validation'])
            mkdir([outputFolderAdjIndx '\Validation']); 
        end
        % To save adjustment index analysis with 5% 
        if ~isdir([outputFolderAdjIndx '\Validation\Ranking_5%'])
            mkdir([outputFolderAdjIndx '\Validation\Ranking_5%']); 
        end
        % To save adjustment index analysis with 10%
        if ~isdir([outputFolderAdjIndx '\Validation\Ranking_10%'])
            mkdir([outputFolderAdjIndx '\Validation\Ranking_10%']); 
        end
        % To save adjustment index analysis with 20%
        if ~isdir([outputFolderAdjIndx '\Validation\Ranking_20%'])
            mkdir([outputFolderAdjIndx '\Validation\Ranking_20%']); 
        end
        % To save adjustment index analysis with 30%          
        if ~isdir([outputFolderAdjIndx '\Validation\Ranking_30%'])
            mkdir([outputFolderAdjIndx '\Validation\Ranking_30%']); 
        end
        % If final Results folder does not exist
        if ~isdir(outputFolderResults)                                   
            mkdir(outputFolderResults);
        end
    end
end
%==========================================================================
% Pre definitions before loop over input prameter lists
%--------------------------------------------------------------------------
% loop counter
counter=0; 
%--------------------------------------------------------------------------
total=ifinal*jfinal*kfinal*lfinal; % total loop length for calculation of 
%                                    run progress (in percentage)
%--------------------------------------------------------------------------
% Creating lists to save outputs
areaList=zeros(ifinal*jfinal*kfinal*lfinal,2); % list for validation via % area
aiList=zeros(ifinal*jfinal*kfinal*lfinal,4); % list for validation via adjust. index
paramsList=zeros(ifinal*jfinal*kfinal*lfinal,4); % list for parameter combinations
czList=[]; % list to check if c/z fraction was already simulated (to avoid 
%            re-work; e.g., c=1000 and z=2 yields same result as c=500 and z=.5)
%==========================================================================
% Checking if program starts from zero
%--------------------------------------------------------------------------
% If so, run proceeds; if no, .bin files are read in order to
% run from where it stopped.
clc; % Clearing last messages
if ~RunInterrupted % If starts from zero
    % Printing message
    fprintf('Done! Initiating calculation...'); 
    fprintf('\n');
    pause(1.5);
    % End of message
else % If run was interrupted (retrieving inputs from .bin files)
    % Printing message
    fprintf('Skiping previous calculations...'); 
    fprintf('\n');
    fprintf('Progress: %s', ['|' [repmat('~',1,round(counter*50/LastCounter+0.1)) '>'] ... 
    	blanks(round((counter-LastCounter)*(50/LastCounter)-0.1)) '| ' num2str(round(100*counter/LastCounter,3),3) '%']);
    fprintf('\n');
    pause(1.5);
    % End of message
    fileBinPA = fopen([outputFolder '\Percent_Area_outputs.bin'],'r');
    areaListOld = fread(fileBinPA,[total 6],'double'); % 6 = 4 params + Min + Mode (Raw not implemented)
    fclose('all');
    %areaListOld = xlsread([outputFolderArea '\Percent_Area_outputs.xlsx'],'sheet1');
    fileBinAI = fopen([outputFolder '\Adjust_Index_outputs.bin'],'r');
    aiListOld = fread(fileBinAI,[total 8],'double'); % 8 = 4 params + 5% + 10% + 20% + 30%
    fclose('all');
    %aiListOld = xlsread([outputFolderAdjIndx '\Adjust_Index_outputs.xlsx'],'sheet1');
    paramsList(1:LastCounter,:)=areaListOld(1:LastCounter,1:4);
    areaList(1:LastCounter,:)=areaListOld(1:LastCounter,5:6);
    aiList(1:LastCounter,:)=aiListOld(1:LastCounter,5:8);
end
%==========================================================================
% Variable to add counter to log file
logInfo=readtable([outputFolder '\log.txt']);
%==========================================================================
% INITIATING LOOP FOR COMBINATION OF INPUT PARAMETERS
%--------------------------------------------------------------------------
% Loop is jumped to last calculation if run was interrupted.
if RunInterrupted && LastCounter == total 
    counter = LastCounter;                
    % Message to indicate that combinations are being skipped
    clc;
    fprintf('Skiping previous calculations...'); 
    fprintf('\n');
    fprintf('Progress: %s', ['|' [repmat('~',1,round(counter*50/LastCounter+0.1)) '>'] ... 
        blanks(round((counter-LastCounter)*(50/LastCounter)-0.1)) '| ' num2str(round(100*counter/LastCounter,3),3) '%']);
    fprintf('\n');
    pause(1.5);
    % End of message
else
    %----------------------------------------------------------------------
    % phi
    for i=1:ifinal 
        if rowsphi==1
            phi = ones(size(THETA1))*phiValues(i)*pi/180; % passing phi to rad
            phiTrack = phiValues(i); % value for paramsList 
        else
            phi = phiValues*pi/180; % passing phi to rad
            phiTrack = NaN; % no value for paramList (case in which input is grid)
        end
        %------------------------------------------------------------------
        % rhoS
        for j=1:jfinal
            if rowsrhoS==1
                rhoS = ones(size(THETA1))*rhoSValues(j);
                rhoSTrack = rhoSValues(j); % value for paramsList 
            else
                rhoS = rhoSvalues;
                rhoSTrack = NaN; % no value for paramList (case in which input is grid)
            end
            %--------------------------------------------------------------
            % c
            for k=1:kfinal
                if rowsc==1
                    c = ones(size(THETA1))*cValues(k);
                    cTrack = cValues(k); % value for paramsList 
                else
                    c = cValues;
                    cTrack = NaN; % no value for paramList (case in which input is grid)
                end
                %----------------------------------------------------------
                % z
                for l=1:lfinal
                    if rowsz==1
                        z = ones(size(THETA1))*zValues(l);
                        zTrack = zValues(l); % value for paramsList 
                    else
                        z = zValues;
                        zTrack = NaN; % no value for paramList (case in which input is grid)
                    end
                    %------------------------------------------------------
                    % Checking if program starts from zero. If was
                    % interrupted, combinations already performed are
                    % skipepd. If starts from zero, calculations proceed.
                    if RunInterrupted == 1 && counter <= (LastCounter -1) 
                        %--------------------------------------------------    
                        % Advancing counter till LastCounter
                        % Printing progreess on screen (only when counter
                        % is divisible by 10)
                        if mod(counter,10)==0
                            clc;
                            fprintf('Skiping previous calculations...');
                            fprintf('\n');
                            fprintf('Progress: %s', ['|' [repmat('~',1,round(counter*50/LastCounter+0.1)) '>'] ... 
                                blanks(round((counter-LastCounter)*(50/LastCounter)-0.1)) '| ' num2str(round(100*counter/LastCounter,3),3) '%']);
                            fprintf('\n');
                        end
                        counter = counter +1;
                        % Delaying final display message (for end of check)
                        if counter == (LastCounter)
                            clc;
                            fprintf('Skiping previous calculations...');
                            fprintf('\n');
                            fprintf('Progress: %s', ['|' [repmat('~',1,round(counter*50/LastCounter+0.1)) '>'] ... 
                                blanks(round((counter-LastCounter)*(50/LastCounter)-0.1)) '| ' num2str(round(100*counter/LastCounter,3),3) '%']);
                            fprintf('\n');
                            pause(1.5);
                        end
                    %------------------------------------------------------
                    % Checking if combination with same c/z was already
                    % calculated.
                    elseif ~isnan(cTrack/zTrack) && ~isempty(czList) && ... % checking if z is list
                            ~isempty(find(czList(:,9)==cTrack/zTrack))     % if czList is empty and c/z is a repeat
                        %--------------------------------------------------
                        % Printing progress 
                        clc; % clearing last messages
                        fprintf('Progress: %s', ['|' [repmat('=',1,round(counter*50/total+0.1)) '>'] ... 
                            blanks(round((total-counter)*(50/total)-0.1)) '| ' num2str(round(100*counter/total,3),3) '%']);
                        fprintf('\n');
                        %--------------------------------------------------
                        % Retrieving position of repeat combination at czList
                        czLoc = find(czList(:,9)==cTrack/zTrack,1);
                        %--------------------------------------------------
                        % Advancing counter
                        counter = counter + 1;
                        %--------------------------------------------------
                        % Savingin counter at log.bin file
                        logBin = fopen([outputFolder '\log.bin'],'w');
                        fileBinColumn(end)=counter;
                        fwrite(logBin, fileBinColumn,'double');
                        fclose('all');
                        % Savingin counter at log.bin file
                        logInfo{14,1}={['Counter ' num2str(counter)]};
                        writetable(logInfo,[outputFolder '\log.txt'],'Delimiter',' ');
                        fclose('all'); 
                        %--------------------------------------------------
                        % Saving parameter combination of current validation
                        paramsList(counter,:) = [czList(czLoc,1) czList(czLoc,2) cTrack zTrack];
                        % Saving validation results for percentage area method
                        areaList(counter,:) = [czList(czLoc,3) czList(czLoc,4)];
                        % Saving validation results for adjustment index method
                        aiList(counter,:) = [czList(czLoc,5) czList(czLoc,6) czList(czLoc,7) czList(czLoc,8)];
                        %--------------------------------------------------
                        % Adding values to Percent_Area_outputs.xlsx file
                        fileBinPA = fopen([outputFolder '\Percent_Area_outputs.bin'],'w');
                        fwrite(fileBinPA, [paramsList areaList],'double');
                        fclose('all');                    
                        % % Adding values to Adjust_Index_outputs.xlsx file
                        fileBinAI = fopen([outputFolder '\Adjust_Index_outputs.bin'],'w');
                        fwrite(fileBinAI, [paramsList aiList],'double');
                        fclose('all');
                        %--------------------------------------------------
                        % Printin progress
                        clc; % clearing last messages
                        fprintf('Progress: %s', ['|' [repmat('=',1,round(counter*50/total+0.1)) '>'] ... 
                            blanks(round((total-counter)*(50/total)-0.1)) '| ' num2str(round(100*counter/total,3),3) '%']);
                        fprintf('\n');    
                    %------------------------------------------------------
                    % If run starts from zero, if combinations are new and
                    % c/z has never being calculated
                    else   
                        % Initiating combination calculations
                        %--------------------------------------------------
                        % Printing progress
                        clc; % clearing last messages
                        fprintf('Progress: %s', ['|' [repmat('=',1,round(counter*50/total+0.1)) '>'] ... 
                            blanks(round((total-counter)*(50/total)-0.1)) '| ' num2str(round(100*counter/total,3),3) '%']);
                        fprintf('\n');
                        %--------------------------------------------------
                        % Calcualting q/T ratio for basin and scars
                        QTBasin = qtgenerator(THETA1, FLOWACC1, phi, rhoS, c, z);
                        QTScars = qtgenerator(THETA2, FLOWACC2, phi, rhoS, c, z);
                        %--------------------------------------------------
                        % Finding scars for uniformization (min and mode).
                        % Calculation is done only on first loop of each
                        % run.
                        if counter==0 || RunInterrupted==1
                            [nscars, ScarsIDs]=scarsidentifier(QTScars);
                            % Finding indexes 
                            index=find(QTScars~=-9999);                        
                            % To avoid calculating again (respecting last if condition)
                            RunInterrupted=0; 
                        end
                        %--------------------------------------------------
                        % Uniformization of q/T in each sacar polygon
                        UnifQTScarsMin = QTScars;  % min
                        UnifQTScarsMode = QTScars; % mode
                        % Finding pixel values inside scars (avoinding NoData)
                        pixels=QTScars(QTScars~=-9999);
                        % Finding list of IDs of each pixel directly from ScarsIDs
                        ID=pixels;
                        for n=1:length(pixels)
                            ID(n)=ScarsIDs(index(n));
                        end
                        % Uniformization
                        for n=1:nscars
                            % Atributing mininum value of a scar to its entirety
                            UnifQTScarsMin(index(ID==n))=min(pixels(ID==n));
                            % % Atributing most repeating value of a scar to its entirety
                            UnifQTScarsMode(index(ID==n))=mode(pixels(ID==n));
                        end
                        %--------------------------------------------------
                        % Calculating integral for percent. area validation  
                        [areaMin,~,~,~,~] = integral(QTBasin, UnifQTScarsMin, resolution);
                        [areaMode,~,~,~,~] = integral(QTBasin, UnifQTScarsMode, resolution);
                        %[areaRaw,~,~,~,~] = integral(QTBasin, QTScars, resolution);
                        %--------------------------------------------------
                        % Avoiding integral errors from some model results
                        % (fail of curve fit or bad performace due to lack
                        % of pixels pertaining to certain classes, creating
                        % gaps)
                        if isnan(areaMin)
                            areaMin=0;                    
                        end
                        %
                        if isnan(areaMode)
                            areaMode=0;
                        end
                        %
                        %if isnan(areaRaw) (not implemented)
                        %    areaRaw=0;
                        %end
                        %--------------------------------------------------
                        % Calculating adjustment index for validation
                        % (ai=adjustment index). This method works only
                        % with q/T grids.
                        ai = adjustindex(QTBasin, QTScars);
                        %--------------------------------------------------
                        % Advancing counter
                        counter = counter + 1;
                        %--------------------------------------------------
                        % Saving counter on log.bin file
                        logBin = fopen([outputFolder '\log.bin'],'w');
                        fileBinColumn(end)=counter;
                        fwrite(logBin, fileBinColumn,'double');
                        fclose('all');
                        % Saving counter on log.txt file
                        logInfo{14,1}={['Counter ' num2str(counter)]};
                        writetable(logInfo,[outputFolder '\log.txt'],'Delimiter',' ');
                        fclose('all'); 
                        %--------------------------------------------------    
                        % Saving parameter combination of current validation
                        paramsList(counter,:) = [phiTrack rhoSTrack cTrack zTrack];
                        % Saving validation values for percent. area method
                        areaList(counter,:) = [areaMin areaMode];
                        % Saving validation values for adjust. index method
                        aiList(counter,:) = ai;                    
                        %--------------------------------------------------                  
                        % Adding values to Percent_Area_outputs.bin file
                        fileBinPA = fopen([outputFolder '\Percent_Area_outputs.bin'],'w');
                        fwrite(fileBinPA, [paramsList areaList],'double');
                        fclose('all');                    
                        % Adding values to Adjust_Index_outputs.bin file
                        fileBinAI = fopen([outputFolder '\Adjust_Index_outputs.bin'],'w');
                        fwrite(fileBinAI, [paramsList aiList],'double');
                        fclose('all');
                        %--------------------------------------------------
                        % Adding c/z value to czList to avoid repeated
                        % calculation (inly if both are inputed as lists)
                        if ~isnan(cTrack) && ~isnan(zTrack)
                            czList(end+1,:) = [phiTrack rhoSTrack areaMin areaMode ai cTrack/zTrack];
                        end
                        %--------------------------------------------------
                        % Printin progress
                        clc; % clearing last messages
                        fprintf('Progress: %s', ['|' [repmat('=',1,round(counter*50/total+0.1)) '>'] ... 
                            blanks(round((total-counter)*(50/total)-0.1)) '| ' num2str(round(100*counter/total,3),3) '%']);
                        fprintf('\n');                    
                    end
                    %------------------------------------------------------
                end % closing z loop
            end % closing c loop
            czList=[]; % rebooting list to track c/z calculation
        end % closing rhoS loop
    end % closing phi loop
    %----------------------------------------------------------------------
    % Saving counter as -9999 to signal that all combinations were
    % calculated and run is complete.
    logBin = fopen([outputFolder '\log.bin'],'w');
    fileBinColumn(end+1)=-9999;
    fwrite(logBin, fileBinColumn,'double');
    fclose('all');
%--------------------------------------------------------------------------
end % closing if before first loop
%==========================================================================
% Message 
clc;
fprintf('Finished! %s', ['|' [repmat('=',1,50) '='] '| ' num2str(100) '%']);
fprintf('\n');                    
fprintf('Saving outputs... 0%%');
fprintf('\n');
pause(1.5);
printercounter=1;
%==========================================================================
% Rounding values of areaList for ordering (similar integral values wont
% allow distinguishing one model's performance from the other)
areaList=round(areaList,6);
%==========================================================================
% Creating .xlsx tables
%--------------------------------------------------------------------------
% Percentage area
xlswrite([outputFolderArea '\Percent_Area_outputs.xlsx'], ...
    {'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '%Area(Min)', '%Area(Mode)'}, 'sheet1', 'A1');
xlswrite([outputFolderArea '\Percent_Area_outputs.xlsx'], ...
    [paramsList areaList], 'sheet1', 'A2');
% Adjustment index
xlswrite([outputFolderAdjIndx '\Adjust_Index_outputs.xlsx'], ...
    {'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '5%', '10%', '20%', '30%'}, 'sheet1', 'A1');
xlswrite([outputFolderAdjIndx '\Adjust_Index_outputs.xlsx'], ...
    [paramsList aiList], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Creating .txt files (with same infos of .xlsx file)
%--------------------------------------------------------------------------
% Opening .txt file for output of percentage area
fclose('all');
fid=fopen([outputFolderArea '\Percent_Area_outputs.txt'], 'w'); 
fprintf(fid, '%s', ['phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '%Area(Min) ' '%Area(Mode) ']);
fprintf(fid, '\n');
for i=1:counter
    fprintf(fid, '%g ', [paramsList(i,:) areaList(i,:)]);
    fprintf(fid, '\n');
end
fclose('all');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing 
%--------------------------------------------------------------------------
% Opening .txt file for output of adjustment index
fclose('all');
fid=fopen([outputFolderAdjIndx '\Adjust_Index_outputs.txt'], 'w'); 
fprintf(fid, '%s', ['phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '5% ' '10% ' '20% ' '30%']);
fprintf(fid, '\n');
for i=1:counter
    fprintf(fid, '%g ', [paramsList(i,:) aiList(i,:)]);
    fprintf(fid, '\n');
end
fclose('all');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Ordering from highest to lowest areaList value
%--------------------------------------------------------------------------
% Min
[areasOrderedMin, indxMin] = sort(areaList(:,1),'descend');
paramsOrderedMin = paramsList(indxMin,:);
% Mode
[areasOrderedMode, indxMode] = sort(areaList(:,2),'descend');
paramsOrderedMode = paramsList(indxMode,:);
% Raw (not implemented)
%[areasOrderedRaw, indxRaw] = sort(areaList(:,3),'descend');
%paramsOrderedRaw = paramsList(indxRaw,:);
% Ordering from higher to lower aiList value for each percentage
% 5
[aiOrdered5, indx5] = sort(aiList(:,1),'descend');
paramsOrdered5 = paramsList(indx5,:);
% 10
[aiOrdered10, indx10] = sort(aiList(:,2),'descend');
paramsOrdered10 = paramsList(indx10,:);
% 20
[aiOrdered20, indx20] = sort(aiList(:,3),'descend');
paramsOrdered20 = paramsList(indx20,:);
% 30
[aiOrdered30, indx30] = sort(aiList(:,4),'descend');
paramsOrdered30 = paramsList(indx30,:);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Ranking best parameter combinations for adjustment index method: the 
% average rank performance of each paramsOrdered ranks (for 5%, 10%, 20% 
% and 30%) is calculated. Next, the averages are ranked from lowest to
% highest since the best performing combination appears closer to top
% positions, e.g., 1 for first, 2 for second, etc.
%--------------------------------------------------------------------------
% List for average of positions 
meanPositions=zeros(size(paramsList,1),1);
posAll=zeros(size(paramsList,1),4);
for i=1:size(aiList,1)
    % 5%
    % Checking if value of a.i. is zero for combination
    if aiList(i,1)==0 % if yes, value of combination position will be same as last
        pos5=size(aiList,1);
    else % if no, value of combination position will be its own position at the ranked i% (i=5, 10..) list 
        pos5=find(indx5==i);
    end
    % 10%
    if aiList(i,2)==0
        pos10=size(aiList,1);
    else
        pos10=find(indx10==i);
    end
    % 20%
    if aiList(i,3)==0
        pos20=size(aiList,1);
    else
        pos20=find(indx20==i);
    end
    % 30%
    if aiList(i,4)==0
        pos30=size(aiList,1);
    else
        pos30=find(indx30==i);
    end
    % Average of combination positions 
    meanPositions(i)=(pos5 + pos10 + pos20 + pos30)/4;
    posAll(i,:) = [pos5 pos10 pos20 pos30];
    %...........................Printing saving progress
    if mod(i,1000)==0 % Printing only when i is divisible by 5
        printercounter=printercounter+1; 
        printer(printercounter,NResults+(2/1000)*size(aiList,1));
    end
    %................................................Printing
end
%--------------------------------------------------------------------------
% Ordering averages from lowest to highest
[Ranks, aiRankIndxs] = sort(meanPositions);
% Creating ordered list (ranked) of position averages
ParamsRanked = paramsList(aiRankIndxs,:);
posAll=posAll(aiRankIndxs,:);
% Creating .txt file with ranks of best performing parameters
fclose('all');
fid=fopen([outputFolderAdjIndx '\Validation\ordered_ai.txt'], 'w'); 
fprintf(fid, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '5% ' '10% ' '20% ' '30% ' 'Mean']);
fprintf(fid, '\n');
%--------------------------------------------------------------------------
% Counter to be used as rank in ordered list
cAI=0;
% List for postiions of combinations in ordered list  
positionsAI=zeros(size(aiList,1),1);
% For checking of repeated parameters
p1AI=0;p2AI=0;p3AI=0;
for i=1:size(aiList,1) % from i=1 to u=n# of lines inside aiList
    % Checking if c/z and other parameters already ranked. If yes, same counter
    % is used to display in ordered list. If no, counter is advanced (+1) for display 
    if sum(([p1AI p2AI p3AI] == [ParamsRanked(i,1:2) ParamsRanked(i,3)/ParamsRanked(i,4)]))~=3
        cAI=cAI+1;
    end
    fprintf(fid, '%g', cAI);
    fprintf(fid, '%s', '. ');
    fprintf(fid, '%g ', [ParamsRanked(i,:) posAll(i,:) Ranks(i,:)]);    
    fprintf(fid, '\n');
    p1AI=ParamsRanked(i,1);
    p2AI=ParamsRanked(i,2);
    p3AI=ParamsRanked(i,3)/ParamsRanked(i,4);
    positionsAI(i)=cAI;
end
fclose('all');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%--------------------------------------------------------------------------
% Creating .xlsx files with rank for best parameters
xlswrite([outputFolderAdjIndx '\Validation\ordered_ai.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]', 'z[m]', '5%', '10%', '20%', '30%', 'Mean'}, 'sheet1', 'A1');
xlswrite([outputFolderAdjIndx '\Validation\ordered_ai.xlsx'], ...
    [positionsAI ParamsRanked posAll Ranks], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Retrieving best performing input parameters for the percentage area method
% (top of lists paramsOrdered for Adjustment index, Min, Mode) (Raw not implemented)
%--------------------------------------------------------------------------
BestParamsAdjIndx = ParamsRanked(1,:);   % Adjust. index
BestParamsMin = paramsOrderedMin(1,:);   % Min
BestParamsMode = paramsOrderedMode(1,:); % Mode
%BestParamsRaw = paramsOrderedRaw(1,:);   % Raw
%==========================================================================
% Ordering lists (area and a.i.) for display in .txt and .xlsx files
%--------------------------------------------------------------------------
% Initiating .txt files
fclose('all');
% Min
fidMin=fopen([outputFolderArea '\Scars_Min\Validation\ordered_Min.txt'], 'w');
fprintf(fidMin, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '%Area(Min)']);
fprintf(fidMin, '\n');
% Mode
fidMode=fopen([outputFolderArea '\Scars_Mode\Validation\ordered_Mode.txt'], 'w');
fprintf(fidMode, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '%Area(Mode)']);
fprintf(fidMode, '\n');
% Raw (not implemented)
%fidRaw=fopen([outputFolderArea '\Scars_Raw\Validation\ordered_Raw.txt'], 'w');
%fprintf(fidRaw, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '%Area(Raw)']);
%fprintf(fidRaw, '\n');
% 5%
fid5=fopen([outputFolderAdjIndx '\Validation\Ranking_5%\ordered_5.txt'], 'w');
fprintf(fid5, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '5%']);
fprintf(fid5, '\n');
% 10%
fid10=fopen([outputFolderAdjIndx '\Validation\Ranking_10%\ordered_10.txt'], 'w');
fprintf(fid10, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '10%']);
fprintf(fid10, '\n');
% 20%
fid20=fopen([outputFolderAdjIndx '\Validation\Ranking_20%\ordered_20.txt'], 'w');
fprintf(fid20, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '20%']);
fprintf(fid20, '\n');
% 30%
fid30=fopen([outputFolderAdjIndx '\Validation\Ranking_30%\ordered_30.txt'], 'w');
fprintf(fid30, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '30%']);
fprintf(fid30, '\n');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%--------------------------------------------------------------------------
% Counters to be used as rank in ordered list
cMin=0; cMode=0; c5=0; c10=0; c20=0; c30=0; % (Raw not implemented)
% Lists for positions of combinations in ordered list
positionsMin=zeros(size(areaList,1),1);
positionsMode=zeros(size(areaList,1),1);
%positionsRaw=zeros(size(areaList,1),1);
positions5=zeros(size(areaList,1),1);
positions10=zeros(size(areaList,1),1);
positions20=zeros(size(areaList,1),1);
positions30=zeros(size(areaList,1),1);
% For checking of repeated parameters
p1Min=0;p2Min=0;p3Min=0;
p1Mode=0;p2Mode=0;p3Mode=0;
%p1Raw=0;p2Raw=0;p3Raw=0;
p15=0;p25=0;p35=0;
p110=0;p210=0;p310=0;
p120=0;p220=0;p320=0;
p130=0;p230=0;p330=0;
%--------------------------------------------------------------------------
% Creating variables for comparing % areas. If a particular % of 2
% combinations are the same, their positions in ordered list is the same.
areaBeforeMin = inf;
areaBeforeMode = inf;
%areaBeforeRaw = inf;
for i=1:size(areaList,1) % from i=1 to i=n# of lines in areaList(=aiList)
    % Min
    % Checking if c/z and other parameters already ranked AND if areaMin(i) is different.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display  
    if sum(([p1Min p2Min p3Min] == [paramsOrderedMin(i,1:2) paramsOrderedMin(i,3)/paramsOrderedMin(i,4)]))~=3 ...
            && areasOrderedMin(i) < areaBeforeMin
        cMin=cMin+1;
    end
    fprintf(fidMin, '%g', cMin);
    fprintf(fidMin, '%s', '. ');
    fprintf(fidMin, '%g ', [paramsOrderedMin(i,:) areasOrderedMin(i,:)]);
    fprintf(fidMin, '\n');
    p1Min=paramsOrderedMin(i,1);
    p2Min=paramsOrderedMin(i,2);
    p3Min=paramsOrderedMin(i,3)/paramsOrderedMin(i,4);
    positionsMin(i)=cMin;
    areaBeforeMin = areasOrderedMin(i);
    %----------------------------------------------------------------------
    % Mode
    % Checking if c/z and other parameters already ranked AND if areaMode(i) is different.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display    
    if sum(([p1Mode p2Mode p3Mode] == [paramsOrderedMode(i,1:2) paramsOrderedMode(i,3)/paramsOrderedMode(i,4)]))~=3 ...
            && areasOrderedMode(i) < areaBeforeMode
        cMode=cMode+1;
    end
    fprintf(fidMode, '%g', cMode);
    fprintf(fidMode, '%s', '. ');
    fprintf(fidMode, '%g ', [paramsOrderedMode(i,:) areasOrderedMode(i,:)]);
    fprintf(fidMode, '\n');
    p1Mode=paramsOrderedMode(i,1);
    p2Mode=paramsOrderedMode(i,2);
    p3Mode=paramsOrderedMode(i,3)/paramsOrderedMode(i,4);
    positionsMode(i)=cMode;
    areaBeforeMode = areasOrderedMode(i);
    %----------------------------------------------------------------------
    % Raw (not implemented)
    % Checking if c/z and other parameters already ranked AND if areaRaw(i) is different.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display  
    %if sum(([p1Raw p2Raw p3Raw] == [paramsOrderedRaw(i,1:2) paramsOrderedRaw(i,3)/paramsOrderedRaw(i,4)]))~=3 ...
    %        && areasOrderedRaw(i) < areaBeforeRaw
    %    cRaw=cRaw+1;
    %end
    %fprintf(fidRaw, '%g', cRaw);
    %fprintf(fidRaw, '%s', '. ');
    %fprintf(fidRaw, '%g ', [paramsOrderedRaw(i,:) areasOrderedRaw(i,:)]);
    %fprintf(fidRaw, '\n');
    %p1Raw=paramsOrderedRaw(i,1);
    %p2Raw=paramsOrderedRaw(i,2);
    %p3Raw=paramsOrderedRaw(i,3)/paramsOrderedRaw(i,4);
    %positionsRaw(i)=cRaw;
    %areaBeforeRaw = areasOrderedRaw(i);
    %----------------------------------------------------------------------
    % 5%
    % Checking if c/z and other parameters already ranked.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display  
    if sum(([p15 p25 p35] == [paramsOrdered5(i,1:2) paramsOrdered5(i,3)/paramsOrdered5(i,4)]))~=3
        c5=c5+1;
    end
    fprintf(fid5, '%g', c5);
    fprintf(fid5, '%s', '. ');
    fprintf(fid5, '%g ', [paramsOrdered5(i,:) aiOrdered5(i,:)]);
    fprintf(fid5, '\n');
    p15=paramsOrdered5(i,1);
    p25=paramsOrdered5(i,2);
    p35=paramsOrdered5(i,3)/paramsOrdered5(i,4);
    positions5(i)=c5;
    %----------------------------------------------------------------------
    % 10%
    % Checking if c/z and other parameters already ranked.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display    
    if sum(([p110 p210 p310] == [paramsOrdered10(i,1:2) paramsOrdered10(i,3)/paramsOrdered10(i,4)]))~=3
        c10=c10+1;
    end
    fprintf(fid10, '%g', c10);
    fprintf(fid10, '%s', '. ');
    fprintf(fid10, '%g ', [paramsOrdered10(i,:) aiOrdered10(i,:)]);
    fprintf(fid10, '\n');
    p110=paramsOrdered10(i,1);
    p210=paramsOrdered10(i,2);
    p310=paramsOrdered10(i,3)/paramsOrdered10(i,4);
    positions10(i)=c10;
    %----------------------------------------------------------------------
    % 20%
    % Checking if c/z and other parameters already ranked.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display     
    if sum(([p120 p220 p320] == [paramsOrdered20(i,1:2) paramsOrdered20(i,3)/paramsOrdered20(i,4)]))~=3
        c20=c20+1;
    end
    fprintf(fid20, '%g', c20);
    fprintf(fid20, '%s', '. ');
    fprintf(fid20, '%g ', [paramsOrdered20(i,:) aiOrdered20(i,:)]);
    fprintf(fid20, '\n');
    p120=paramsOrdered20(i,1);
    p220=paramsOrdered20(i,2);
    p320=paramsOrdered20(i,3)/paramsOrdered20(i,4);
    positions20(i)=c20;
    %----------------------------------------------------------------------
    % 30%
    % Checking if c/z and other parameters already ranked.
    % If yes, same counter is used to display in ordered list. If no, 
    % counter is advanced (+1) for display     
    if sum(([p130 p230 p330] == [paramsOrdered30(i,1:2) paramsOrdered30(i,3)/paramsOrdered30(i,4)]))~=3
        c30=c30+1;
    end
    fprintf(fid30, '%g', c30);
    fprintf(fid30, '%s', '. ');
    fprintf(fid30, '%g ', [paramsOrdered30(i,:) aiOrdered30(i,:)]);
    fprintf(fid30, '\n');
    p130=paramsOrdered30(i,1);
    p230=paramsOrdered30(i,2);
    p330=paramsOrdered30(i,3)/paramsOrdered30(i,4);
    positions30(i)=c30;
    %...........................Printing saving progress
    if mod(i,1000)==0 % Printing only when i is divisible by 5
        printercounter=printercounter+1; 
        printer(printercounter,NResults+(2/1000)*size(aiList,1));
    end
    %................................................Printing
end
fclose('all');
%--------------------------------------------------------------------------
% Wiriting .xlsx files
fclose('all');
% Min
xlswrite([outputFolderArea '\Scars_Min\Validation\ordered_Min.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '%Area(Min)'}, 'sheet1', 'A1');
xlswrite([outputFolderArea '\Scars_Min\Validation\ordered_Min.xlsx'], ...
    [positionsMin paramsOrderedMin areasOrderedMin], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Mode
xlswrite([outputFolderArea '\Scars_Mode\Validation\ordered_Mode.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '%Area(Mode)'}, 'sheet1', 'A1');
xlswrite([outputFolderArea '\Scars_Mode\Validation\ordered_Mode.xlsx'], ...
    [positionsMode paramsOrderedMode areasOrderedMode], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Raw (not implemented)
%xlswrite([outputFolderArea '\Scars_Raw\Validation\ordered_Raw.xlsx'], ...
%    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '%Area(Raw)'}, 'sheet1', 'A1');
%xlswrite([outputFolderArea '\Scars_Raw\Validation\ordered_Raw.xlsx'], ...
%    [positionsRaw paramsOrderedRaw areasOrderedRaw], 'sheet1', 'A2');
%...........................Printing saving progress
%printercounter=printercounter+1; 
%printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% 5%
xlswrite([outputFolderAdjIndx '\Validation\Ranking_5%\ordered_5.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '5%'}, 'sheet1', 'A1');
xlswrite([outputFolderAdjIndx '\Validation\Ranking_5%\ordered_5.xlsx'], ...
    [positions5 paramsOrdered5 aiOrdered5], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% 10%
xlswrite([outputFolderAdjIndx '\Validation\Ranking_10%\ordered_10.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '10%'}, 'sheet1', 'A1');
xlswrite([outputFolderAdjIndx '\Validation\Ranking_10%\ordered_10.xlsx'], ...
    [positions10 paramsOrdered10 aiOrdered10], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% 20%
xlswrite([outputFolderAdjIndx '\Validation\Ranking_20%\ordered_20.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '20%'}, 'sheet1', 'A1');
xlswrite([outputFolderAdjIndx '\Validation\Ranking_20%\ordered_20.xlsx'], ...
    [positions20 paramsOrdered20 aiOrdered20], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% 30%
xlswrite([outputFolderAdjIndx '\Validation\Ranking_30%\ordered_30.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '30%'}, 'sheet1', 'A1');
xlswrite([outputFolderAdjIndx '\Validation\Ranking_30%\ordered_30.xlsx'], ...
    [positions30 paramsOrdered30 aiOrdered30], 'sheet1', 'A2');
fclose('all');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Working basin and scars susceptibility maps for figures
%--------------------------------------------------------------------------
% Creating data cube of 2D stacked data (n rows, ncols, 2D grid of input
% value
[rows, cols] = size(THETA1);
%DataCubeParamsRaw = zeros(rows, cols, 4); (not implemented)
DataCubeParamsMin = zeros(rows, cols, 4);
DataCubeParamsMode = zeros(rows, cols, 4);
DataCubeParamsAdjIndx = zeros(rows, cols, 4);
%--------------------------------------------------------------------------
% Input list to be used on loop
Args = {phiValues rhoSValues cValues zValues};
for i=1:4
    if isnan(BestParamsMin(i))        
        DataCubeParamsMin(:,:,i)=Args{i};  % If input is grid (Min)
        DataCubeParamsMode(:,:,i)=Args{i}; %         ''       (Mode)
        %DataCubeParamsRaw(:,:,i)=Args{i};  %         ''       (Raw)
        DataCubeParamsAdjIndx(:,:,i)=Args{i};  %     ''       (AdjIndx) 
    else        
        DataCubeParamsMin(:,:,i)=ones(size(THETA1))*BestParamsMin(i);    %  If list: (Min)
        DataCubeParamsMode(:,:,i)=ones(size(THETA1))*BestParamsMode(i);  %      ''   (Mode)
        %DataCubeParamsRaw(:,:,i)=ones(size(THETA1))*BestParamsRaw(i);    %      ''   (Raw)
        DataCubeParamsAdjIndx(:,:,i)=ones(size(THETA1))*BestParamsAdjIndx(i); % ''   (AdjIndx)
    end
end
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing 
%==========================================================================
% Using grids of data cube as argument inside figure making functions;
% first best maps of q/T and shalstab classes for both basin and scars is
% selected
%--------------------------------------------------------------------------
% Min 
BestQTBasinMin = qtgenerator(THETA1, FLOWACC1, ... 
    DataCubeParamsMin(:,:,1)*pi/180, DataCubeParamsMin(:,:,2), DataCubeParamsMin(:,:,3), DataCubeParamsMin(:,:,4));
BestQTScarsMin = qtgenerator(THETA2, FLOWACC2, ... 
    DataCubeParamsMin(:,:,1)*pi/180, DataCubeParamsMin(:,:,2), DataCubeParamsMin(:,:,3), DataCubeParamsMin(:,:,4));
BestClassBasinMin = classgenerator(THETA1, FLOWACC1, ... 
    DataCubeParamsMin(:,:,1)*pi/180, DataCubeParamsMin(:,:,2), DataCubeParamsMin(:,:,3), DataCubeParamsMin(:,:,4));
BestClassScarsMin = classgenerator(THETA2, FLOWACC2, ... 
    DataCubeParamsMin(:,:,1)*pi/180, DataCubeParamsMin(:,:,2), DataCubeParamsMin(:,:,3), DataCubeParamsMin(:,:,4));
% Making scar maps uniform (Min)
[nscars, ScarsIDs]=scarsidentifier(BestQTScarsMin);
UnifQTScarsMin = BestQTScarsMin;
%UnifClassScarsMin = BestClassScarsMin;
for id=1:nscars
    UnifQTScarsMin(ScarsIDs==id)=min(BestQTScarsMin(ScarsIDs==id));
    %UnifClassScarsMin(ScarsIDs==id)=min(BestClassScarsMin(ScarsIDs==id));
end
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%--------------------------------------------------------------------------
% Mode
BestQTBasinMode = qtgenerator(THETA1, FLOWACC1, ... 
    DataCubeParamsMode(:,:,1)*pi/180, DataCubeParamsMode(:,:,2), DataCubeParamsMode(:,:,3), DataCubeParamsMode(:,:,4));
BestQTScarsMode = qtgenerator(THETA2, FLOWACC2, ... 
    DataCubeParamsMode(:,:,1)*pi/180, DataCubeParamsMode(:,:,2), DataCubeParamsMode(:,:,3), DataCubeParamsMode(:,:,4));
BestClassBasinMode = classgenerator(THETA1, FLOWACC1, ... 
    DataCubeParamsMode(:,:,1)*pi/180, DataCubeParamsMode(:,:,2), DataCubeParamsMode(:,:,3), DataCubeParamsMode(:,:,4));
BestClassScarsMode = classgenerator(THETA2, FLOWACC2, ... 
    DataCubeParamsMode(:,:,1)*pi/180, DataCubeParamsMode(:,:,2), DataCubeParamsMode(:,:,3), DataCubeParamsMode(:,:,4));
% Making scar maps uniform (Mode)
[nscars, ScarsIDs]=scarsidentifier(BestQTScarsMode);
UnifQTScarsMode = BestQTScarsMode;
%UnifClassScarsMode = BestClassScarsMode;
for id=1:nscars
    UnifQTScarsMode(ScarsIDs==id)=mode(BestQTScarsMode(ScarsIDs==id));
    %UnifClassScarsMode(ScarsIDs==id)=mode(BestClassScarsMode(ScarsIDs==id));
end
%...................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%--------------------------------------------------------------------------
% Raw (not implemented)
%BestQTBasinRaw = qtgenerator(THETA1, FLOWACC1, ... 
%    DataCubeParamsRaw(:,:,1)*pi/180, DataCubeParamsRaw(:,:,2), DataCubeParamsRaw(:,:,3), DataCubeParamsRaw(:,:,4));
%BestQTScarsRaw = qtgenerator(THETA2, FLOWACC2, ... 
%    DataCubeParamsRaw(:,:,1)*pi/180, DataCubeParamsRaw(:,:,2), DataCubeParamsRaw(:,:,3), DataCubeParamsRaw(:,:,4));
%BestClassBasinRaw = classgenerator(THETA1, FLOWACC1, ... 
%    DataCubeParamsRaw(:,:,1)*pi/180, DataCubeParamsRaw(:,:,2), DataCubeParamsRaw(:,:,3), DataCubeParamsRaw(:,:,4));
%BestClassScarsRaw = classgenerator(THETA2, FLOWACC2, ... 
%    DataCubeParamsRaw(:,:,1)*pi/180, DataCubeParamsRaw(:,:,2), DataCubeParamsRaw(:,:,3), DataCubeParamsRaw(:,:,4));
%...........................Printing saving progress
%printercounter=printercounter+1; 
%printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%--------------------------------------------------------------------------
% AdjIndx
BestQTBasinAdjIndx = qtgenerator(THETA1, FLOWACC1, ... 
    DataCubeParamsAdjIndx(:,:,1)*pi/180, DataCubeParamsAdjIndx(:,:,2), DataCubeParamsAdjIndx(:,:,3), DataCubeParamsAdjIndx(:,:,4));
BestQTScarsAdjIndx = qtgenerator(THETA2, FLOWACC2, ... 
    DataCubeParamsAdjIndx(:,:,1)*pi/180, DataCubeParamsAdjIndx(:,:,2), DataCubeParamsAdjIndx(:,:,3), DataCubeParamsAdjIndx(:,:,4));
BestClassBasinAdjIndx = classgenerator(THETA1, FLOWACC1, ... 
    DataCubeParamsAdjIndx(:,:,1)*pi/180, DataCubeParamsAdjIndx(:,:,2), DataCubeParamsAdjIndx(:,:,3), DataCubeParamsAdjIndx(:,:,4));
BestClassScarsAdjIndx = classgenerator(THETA2, FLOWACC2, ... 
    DataCubeParamsAdjIndx(:,:,1)*pi/180, DataCubeParamsAdjIndx(:,:,2), DataCubeParamsAdjIndx(:,:,3), DataCubeParamsAdjIndx(:,:,4));
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Creating the names of output files
%--------------------------------------------------------------------------
% Min
fileNameMin = ['_phi_' num2str(BestParamsMin(1)) '_rhoS_' num2str(BestParamsMin(2)) ...
    '_c_' num2str(BestParamsMin(3)) '_z_' num2str(BestParamsMin(4))];
% Mode
fileNameMode = ['_phi_' num2str(BestParamsMode(1)) '_rhoS_' num2str(BestParamsMode(2)) ...
    '_c_' num2str(BestParamsMode(3)) '_z_' num2str(BestParamsMode(4))];
% Raw (not implemented)
%fileNameRaw = ['_phi_' num2str(BestParamsRaw(1)) '_rhoS_' num2str(BestParamsRaw(2)) ...
%    '_c_' num2str(BestParamsRaw(3)) '_z_' num2str(BestParamsRaw(4))];
% AdjIndx
fileNameAdjIndx = ['_phi_' num2str(BestParamsAdjIndx(1)) '_rhoS_' num2str(BestParamsAdjIndx(2)) ...
    '_c_' num2str(BestParamsAdjIndx(3)) '_z_' num2str(BestParamsAdjIndx(4))];
%==========================================================================
% Creating ascii grid files of the suceptibility maps for both basin and scars
%--------------------------------------------------------------------------
% Min
asciimaker(thetamap1, BestQTBasinMin, [outputFolderArea '\Scars_Min\Grids'], ['qt_ratio_basin_min_' fileNameMin]);
asciimaker(thetamap2, BestQTScarsMin, [outputFolderArea '\Scars_Min\Grids'], ['qt_ratio_scars_min_' fileNameMin]);
asciimaker(thetamap1, BestClassBasinMin, [outputFolderArea '\Scars_Min\Grids'], ['suscept_map_basin_min_' fileNameMin]);
asciimaker(thetamap2, BestClassScarsMin, [outputFolderArea '\Scars_Min\Grids'], ['suscept_map_scars_min_' fileNameMin]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Mode
asciimaker(thetamap1, BestQTBasinMode, [outputFolderArea '\Scars_Mode\Grids'], ['qt_ratio_basin_mode_' fileNameMode]);
asciimaker(thetamap2, BestQTScarsMode, [outputFolderArea '\Scars_Mode\Grids'], ['qt_ratio_map_scars_mode_' fileNameMode]);
asciimaker(thetamap1, BestClassBasinMode, [outputFolderArea '\Scars_Mode\Grids'], ['suscept_map_basin_mode_' fileNameMode]);
asciimaker(thetamap2, BestClassScarsMode, [outputFolderArea '\Scars_Mode\Grids'], ['suscept_map_scars_mode_' fileNameMode]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Raw (not implemented)
%asciimaker(thetamap1, BestQTBasinRaw, [outputFolderArea '\Scars_Raw\Grids'], ['qt_ratio_basin_raw_' fileNameRaw]);
%asciimaker(thetamap2, BestQTScarsRaw, [outputFolderArea '\Scars_Raw\Grids'], ['qt_ratio_scars_raw_' fileNameRaw]);
%asciimaker(thetamap1, BestClassBasinRaw, [outputFolderArea '\Scars_Raw\Grids'], ['suscept_map_basin_raw_' fileNameRaw]);
%asciimaker(thetamap2, BestClassScarsRaw, [outputFolderArea '\Scars_Raw\Grids'], ['suscept_map_scars_raw_' fileNameRaw]);
%...........................Printing saving progress
%printercounter=printercounter+1; 
%printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% AdjIndx
asciimaker(thetamap1, BestQTBasinAdjIndx, [outputFolderAdjIndx '\Grids'], ['qt_ratio_basin_ai_' fileNameAdjIndx]);
asciimaker(thetamap2, BestQTScarsAdjIndx, [outputFolderAdjIndx '\Grids'], ['qt_ratio_scars_ai_' fileNameAdjIndx]);
asciimaker(thetamap1, BestClassBasinAdjIndx, [outputFolderAdjIndx '\Grids'], ['suscept_map_basin_ai_' fileNameAdjIndx]);
asciimaker(thetamap2, BestClassScarsAdjIndx, [outputFolderAdjIndx '\Grids'], ['suscept_map_scars_ai_' fileNameAdjIndx]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Creating class distribution plots, barplot and integral plots
%--------------------------------------------------------------------------
% Min
classdistplotmaker(THETA1, FLOWACC1, BestClassBasinMin, ...
    [outputFolderArea '\Scars_Min\Sucept_Plots'], ['min_' fileNameMin]);
barplotmaker(BestClassBasinMin, BestClassScarsMin, ... 
    [outputFolderArea '\Scars_Min\Sucept_Plots'], ['min_' fileNameMin]);
areaplotmaker(BestQTBasinMin, UnifQTScarsMin, resolution, ... 
    [outputFolderArea '\Scars_Min\Validation'], ['min_' fileNameMin]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Mode
classdistplotmaker(THETA1, FLOWACC1, BestClassBasinMode, ...
    [outputFolderArea '\Scars_Mode\Sucept_Plots'], ['mode_' fileNameMode]);
barplotmaker(BestClassBasinMode, BestClassScarsMode, ... 
    [outputFolderArea '\Scars_Mode\Sucept_Plots'], ['mode_' fileNameMode]);
areaplotmaker(BestQTBasinMode, UnifQTScarsMode, resolution, ... 
    [outputFolderArea '\Scars_Mode\Validation'], ['mode_' fileNameMode]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Raw (not implemented)
%classdistplotmaker(THETA1, FLOWACC1, BestClassBasinRaw, ...
%    [outputFolderArea '\Scars_Raw\Sucept_Plots'], ['raw_' fileNameRaw]);
%barplotmaker(BestClassBasinRaw, BestClassScarsRaw, ... 
%    [outputFolderArea '\Scars_Raw\Sucept_Plots'], ['raw_' fileNameRaw]);
%areaplotmaker(BestQTBasinRaw, BestQTScarsRaw, resolution, ... 
%    [outputFolderArea '\Scars_Raw\Validation'], ['raw_' fileNameRaw]);
%...........................Printing saving progress
%printercounter=printercounter+1;
%printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% AdjIndx
classdistplotmaker(THETA1, FLOWACC1, BestClassBasinAdjIndx, ...
    [outputFolderAdjIndx '\Sucept_Plots'], ['ai_' fileNameAdjIndx]);
barplotmaker(BestClassBasinAdjIndx, BestClassScarsAdjIndx, ... 
    [outputFolderAdjIndx '\Sucept_Plots'], ['ai_' fileNameAdjIndx]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Creating figures of the suceptibility maps for both basin and scars 
%--------------------------------------------------------------------------
% Min
classfigsmaker(BestClassBasinMin, [outputFolderArea '\Scars_Min\Sucept_Plots'], ['suscept_map_basin_min_' fileNameMin]);
classfigsmaker(BestClassScarsMin, [outputFolderArea '\Scars_Min\Sucept_Plots'], ['suscept_map_scars_min_' fileNameMin]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Mode
classfigsmaker(BestClassBasinMode, [outputFolderArea '\Scars_Mode\Sucept_Plots'], ['suscept_map_basin_mode_' fileNameMode]);
classfigsmaker(BestClassScarsMode, [outputFolderArea '\Scars_Mode\Sucept_Plots'], ['suscept_map_scars_mode_' fileNameMode]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% Raw (not implemented)
%classfigsmaker(BestClassBasinRaw, [outputFolderArea '\Scars_Raw\Sucept_Plots'], ['suscept_map_basin_raw_' fileNameRaw]);
%classfigsmaker(BestClassScarsRaw, [outputFolderArea '\Scars_Raw\Sucept_Plots'], ['suscept_map_scars_raw_' fileNameRaw]);
%...........................Printing saving progress
%printercounter=printercounter+1; 
%printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
% AdjIndx
classfigsmaker(BestClassBasinAdjIndx, [outputFolderAdjIndx '\Sucept_Plots'], ['suscept_map_basin_ai_' fileNameAdjIndx]);
classfigsmaker(BestClassScarsAdjIndx, [outputFolderAdjIndx '\Sucept_Plots'], ['suscept_map_scars_ai_' fileNameAdjIndx]);
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Creating final result list of each param. combination performance inside 
% orderedList (Min, Mode and A.I.) - this is the general performance rank.
%--------------------------------------------------------------------------
FinalOrder=zeros(size(paramsList,1), 3);
for i=1:size(areaList,1) % from i=1 to i=n# of rows in areaList(=aiList)
    % Min
    FinalOrder(indxMin(i),1)=positionsMin(i);
    % Mode
    FinalOrder(indxMode(i),2)=positionsMode(i);
    % AI
    FinalOrder(aiRankIndxs(i),3)=positionsAI(i);    
end
% Ordering final list from lowest to highest with average of the 3 cols
[OrderMeanValues, OrderMeanRank]=sort(mean(FinalOrder,2));
% Ordering list of values for Min, Mode and AI from lowest orderRank to
% highest
AllOrders=FinalOrder(OrderMeanRank,:);
% Ordering parameter list from lowest orderRank to highest
paramsOrderFinal = paramsList(OrderMeanRank,:);
% If ranking of two param. combinations are the same, their position in the
% final list is the same
cOrder = 0;
OrderMeanValueBefore = -inf;
% List for position of combinations in final ordered list
positionsOrder=zeros(size(areaList,1),1);
% Initiating final list of ranks
fclose('all');
fid=fopen([outputFolderResults '\final_result.txt'], 'w'); 
fprintf(fid, '%s', ['Pos ' 'phi[deg] ' 'rhoS[kg/m3] ' 'c[Pa] ' 'z[m] ' '%Area(Min) ' '%Area(Mode) ' 'AI ' 'Mean']);
fprintf(fid, '\n');
for i=1:counter
    if OrderMeanValues(i) > OrderMeanValueBefore
        cOrder=cOrder+1;
    end
    fprintf(fid, '%g', cOrder);
    fprintf(fid, '%s', '. ');
    fprintf(fid, '%g ', [paramsOrderFinal(i,:) AllOrders(i,:) OrderMeanValues(i,:)]);
    fprintf(fid, '\n');
    positionsOrder(i) = cOrder;
    OrderMeanValueBefore = OrderMeanValues(i);
end
OrderMeanValues(i)
fclose('all'); 
% Writing the same resulting list in a final .xlsx file
fclose('all');
xlswrite([outputFolderResults '\final_result.xlsx'], ...
    {'Pos', 'phi[deg]', 'rhoS[kg/m3]', 'c[Pa]' 'z[m]', '%Area(Min)', '%Area(Mode)', 'AI', 'Mean'}, 'sheet1', 'A1');
xlswrite([outputFolderResults '\final_result.xlsx'], ...
    [positionsOrder paramsOrderFinal AllOrders OrderMeanValues], 'sheet1', 'A2');
%...........................Printing saving progress
printercounter=printercounter+1; 
printer(printercounter,NResults+(2/1000)*size(aiList,1));
%................................................Printing
%==========================================================================
% Creating outputs of best performing parameter combinations
%--------------------------------------------------------------------------
% Creating data cube of 2D stacked data (n rows, ncols, 2D grid of input
% value)
[rows, cols] = size(THETA1);
DataCubeParams = zeros(rows, cols, 4);
%Input list to be used on loop
Args = {phiValues rhoSValues cValues zValues};
% Initiating loop
i=1; % counter
% Working NResults, number of models to be considered for output generation
NRLoop = NResults;
if NRLoop > counter % If number of results chosen by user is > n# of possible 
    NRLoop = 1;     % combinaionts (by mistake), NResults is automatically set to 1.
end                   
check=0; % to stop while loop
previousOMV = 0; % to check if model value is different from last (avoiding
                 % outputs with same c/z.)
while check<NRLoop
    if OrderMeanValues(i) ~= previousOMV
        params = [paramsOrderFinal(i,1) paramsOrderFinal(i,2) paramsOrderFinal(i,3) paramsOrderFinal(i,4)];
        % Generating all outputs for top parameter combinations
        fileName = ['phi_' num2str(params(1)) '_rhoS_' num2str(params(2)) '_c_' num2str(params(3)) '_z_' num2str(params(4))];    
        folderName=[outputFolderResults '\' fileName];
        % Checking if output folder exists
        if isdir(folderName)     
            rmdir(folderName, 's');  % If yes, it is deleted and recreated
            mkdir(folderName);
        else
            mkdir(folderName);  % If no, it is created
        end
        % Defyning data cube values with each parameter combination entry
        for j=1:4
            if isnan(params(j))        
                DataCubeParams(:,:,j)=Args{j};  % if input is grid
            else        
                DataCubeParams(:,:,j)=ones(size(THETA1))*params(j);    %  if list
            end
        end
        % Creating maps for the chosen parameter combination 
        BestQTBasin = qtgenerator(THETA1, FLOWACC1, DataCubeParams(:,:,1)*pi/180, ...
            DataCubeParams(:,:,2), DataCubeParams(:,:,3), DataCubeParams(:,:,4));
        BestQTScars = qtgenerator(THETA2, FLOWACC2, DataCubeParams(:,:,1)*pi/180, ...
            DataCubeParams(:,:,2), DataCubeParams(:,:,3), DataCubeParams(:,:,4));
        BestClassBasin = classgenerator(THETA1, FLOWACC1, DataCubeParams(:,:,1)*pi/180, ...
            DataCubeParams(:,:,2), DataCubeParams(:,:,3), DataCubeParams(:,:,4));
        BestClassScars = classgenerator(THETA2, FLOWACC2, DataCubeParams(:,:,1)*pi/180, ...
            DataCubeParams(:,:,2), DataCubeParams(:,:,3), DataCubeParams(:,:,4));
        % Creating ascii grid files of the suceptibility maps for both
        % basin and scars for the chosen parameter combination 
        asciimaker(thetamap1, BestQTBasin, folderName, ['qt_ratio_basin__' fileName]);
        asciimaker(thetamap2, BestQTScars, folderName, ['qt_ratio_scars__' fileName]);
        asciimaker(thetamap1, BestClassBasin, folderName, ['suscept_map_basin__' fileName]);
        asciimaker(thetamap2, BestClassScars, folderName, ['suscept_map_scars__' fileName]);
        % Creating class distribution plots, barplot and integral plots for
        % the chosen parameter combination 
        classdistplotmaker(THETA1, FLOWACC1, BestClassBasin, folderName, ['_' fileName]);
        barplotmaker(BestClassBasin, BestClassScars, folderName, ['_' fileName]);   
        % Creating figures of the suceptibility maps for both basin and 
        % scars for the chosen parameter combination 
        classfigsmaker(BestClassBasin, folderName, ['suscept_map_basin_min__' fileName]);
        classfigsmaker(BestClassScars, folderName, ['suscept_map_scars_min__' fileName]);
        % Advancing loop check
        check = positionsOrder(i);
        %...........................Printing saving progress
        printercounter=printercounter+1; 
        printer(printercounter,NResults+(2/1000)*size(aiList,1));
        %................................................Printing
    end
    % Registering previous average
    previousOMV = OrderMeanValues(i);
    % Advancing counter
    i = i+1;  
end
%==========================================================================
% Wrapping up last lines of lo.bin file (adding finishing time) and
% altering Counter to complete
%--------------------------------------------------------------------------
logInfo{14,1}={['Counter ' 'complete']};
writetable(logInfo,[outputFolder '\log.txt'],'Delimiter',' ');
%
logInfo{16,1}={['Finished at ' datestr(datetime('now'))]};
writetable(logInfo,[outputFolder '\log.txt'],'Delimiter',' ');
%==========================================================================
% Ending script!
%.........................Printing saving progress
printercounterFinal=-9999; 
printer(printercounterFinal,NResults+(2/5)*size(aiList,1));
%..............................................Printing
fclose('all');
clc; % clearing last message
fprintf('All done!');
fprintf('\n');
end