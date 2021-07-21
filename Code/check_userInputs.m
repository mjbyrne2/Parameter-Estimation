%% check_userInputs.m
% Checks the inputs specified for use in AdaptiveWindowedRegularization.m.
% This script can call prompt_userInputs.m if no inputs have been 
% specified.
%
% check_userInputs.m uses penaltyMatrix.m to generate the list of valid
% penalty matrices.

% Default user inputs:
defaultUserInputs.blur = 16;    % Medium blur
defaultUserInputs.SNR = 25;
defaultUserInputs.penalty = "Identity";
defaultUserInputs.windows = {2,"linear"}; % Two linearly spaced windows
defaultUserInputs.resolutions = 0; % Full problem (256 x 256)
defaultUserInputs.methods = ["Best","UPRE","GCV"]; % Parameter selection methods
defaultUserInputs.shuffle = []; % Shuffle images
% <-- Can add more fields if necessary

% Valid options:
validPenalties = penaltyMatrix();   % Get valid penalty matrices
validWindows = ["linear","linearCosine","log","logCosine"];
validMethods = ["Best","UPRE","GCV","MDP"];
validResolutions = [0,1,2,3,4,5];

% Get field names from default inputs:
userInputsFields = string(fieldnames(defaultUserInputs));

inputInstructions = 'Enter ''1'' for default inputs, ''2'' for user specified inputs, or ''0'' to cancel: ';
if ~exist('userInputs','var')
    userResponse = input(['No user inputs have been specified for AdaptiveWindowedRegularization.m.',...
        newline,inputInstructions]);
    switch userResponse
        case 0
            error('AdaptiveWindowedRegularization.m has been halted by the user.')
        case 1
            userInputs = defaultUserInputs;
            disp('The fields of userInputs have been set to default:')
            disp(struct2table(userInputs,'AsArray',true))
        case 2
            run('prompt_userInputs.m')
        otherwise
            error('Invalid response. AdaptiveWindowedRegularization.m has been halted.')
    end
elseif ~isa(userInputs,'struct')    % Check if userInputs is not of type struct
    userResponse = input('The variable userInputs exists but is not of type struct. Do you want to redefine userInputs and continue? (Y/N) \n','s');
    if strcmp(userResponse,'Y')
        clear userInputs
        run('prompt_userInputs.m')
    else
        error('AdaptiveWindowedRegularization.m has been halted by the user.')
    end
elseif ~isempty(setdiff(userInputsFields,string(fieldnames(userInputs))))
    missingFields = setdiff(userInputsFields,string(fieldnames(userInputs)));
    error(strcat("Error: The following fields are missing from userInputs: ",...
        strjoin(missingFields),". The missing fields must be added before AdaptiveWindowedRegularization.m will run."))
else
    disp(strcat("Caution: The fields of userInputs might not have proper syntax.",...
        newline,"See prompt_userInputs.m and check_userInputs.m for information regarding field syntax."))
end

%% Check the certain fields of userInputs are valid for later use

% Check penalty matrix:
if isempty(userInputs.penalty) || ~ismember(userInputs.penalty,validPenalties)
    error('Error: The input penalty matrix is invalid.')
end

% Check windows:
if numel(userInputs.windows) == 2 && ~ismember(userInputs.windows{2},validWindows)
    error(strcat("Error: Invalid window type. Supported window types: ",strjoin(validWindows)))
end

% Check methods:
if ~all(ismember(userInputs.methods,validMethods))
    error(strcat("Error: Invalid parameter selection method(s) detected: ",...
        strjoin(setdiff(userInputs.methods,validMethods)),". Supported methods: ",strjoin(validMethods)))
end

% Check downsampling resolutions:
if ~all(ismember(userInputs.resolutions,validResolutions))
    error(strcat("Error: Invalid downsampling resolution(s) detected: ",...
        strjoin(string(setdiff(userInputs.resolutions,validResolutions))),". Supported resolutions: ",strjoin(string(validResolutions))))
end
