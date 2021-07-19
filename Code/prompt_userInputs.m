%% prompt_userInputs.m
% This script prompts the user for inputs that will be stored as the
% structure "userInputs" and used in AdaptiveWindowedRegularization.m. The
% script is not guaranteed to run if a variable called userInputs is
% already defined; in such case, it is recommended that the variable be
% cleared before running this script.
%
% userInputs has the following fields
% - userInputs.blur (required): Non-negative scalar
% - userInputs.SNR (required): Non-negative scalar or interval.
% - userInputs.penalty (required): String. If the user does not enter a
% string, the default value of userInputs.penalty is empty.
% - userInputs.windows (optional): Cell array of the form 
% {positive integer, non-empty string}. If the user's format does not 
% match, the default value of userInputs.windows is empty.  
% - userInputs.resolutions (optional): Vector of non-negative integers. Any
% duplicate elements are removed. If the user's format does not match, the 
% default value of userInputs.windows is empty.

disp('The blur amount, SNR, and penalty matrix must be specified. The windows and downsampling resolutions can be empty.')

% Initialize empty struture fields:
userInputs.blur = [];
userInputs.SNR = [];
userInputs.penalty = [];
userInputs.windows = [];
userInputs.resolutions = [];
userInputs.methods = [];
userInputs.shuffle = [];

% Set the amount of image blur:
blurClasses = {'numeric'};
blurAttributes = {'nonnegative','scalar','nonempty'};
while isempty(userInputs.blur)
    userInputs.blur = input(['Specify amount of image blur (non-negative number).',...
        newline,'(Enter 9 for low blur, 16 for medium blur, and 25 for severe blur): ']);
    try
        validateattributes(userInputs.blur,blurClasses,blurAttributes)
    catch
        disp('Error: The amount of images blur must be a non-negative number.')
        userInputs.blur = [];
    end 
end
clear blurClasses blurAttributes

% Set the desired SNR of data to be generated:
SNRClasses = {'numeric'};
SNRAttributes = {'nonnegative','vector','nonempty'};
while isempty(userInputs.SNR)
    userInputs.SNR = input(['Specify SNR of data to be generated (non-negative number or interval).',...
        newline,'(Enter 10 for high noise, 25 for medium noise, and 40 for low noise): ']);
    try
        validateattributes(userInputs.SNR,SNRClasses,SNRAttributes)
    catch
        disp('Error: The SNR of the data to be generated must be a non-negative number or interval.')
        userInputs.SNR = [];
    end
    
    % Check size of SNR input and sort if needed:
    if numel(userInputs.SNR) == 1
        break
    elseif numel(userInputs.SNR) > 2
        disp('Error: The SNR of the data to be generated must be a non-negative number or interval.')
        userInputs.SNR = [];     
    elseif numel(userInputs.SNR) == 2
        userInputs.SNR = sort(userInputs.SNR,'ascend');
        if userInputs.SNR(1) == userInputs.SNR(2)
            userInputs.SNR = userInputs.SNR(1);
            disp(['Warning: Interval of measure zero converted to the single SNR of ',...
            num2str(userInputs.SNR),'.'])
        end
    end
    
end
clear SNRClasses SNRAttributes

% Set the penalty matrix:
penaltyClasses = {'string'};
penaltyAttributes = {'scalartext','nonempty'};
while isempty(userInputs.penalty)
    userInputs.penalty = input(['Specify penalty matrix (string).',newline,...
        '(Supported penalty matrices are ',...
        convertStringsToChars(strcat("""",join(validPenalties,""", """),"""")),...
        '): ']);
    try
        validateattributes(userInputs.penalty,penaltyClasses,penaltyAttributes)
    catch
        disp('Error: The penalty matrix must be input as a string.')
        userInputs.penalty = [];
    end 
end
clear penaltyClasses penaltyAttributes

% Set the number and type of spectral windows:
% Example: {2,"logCosine"}
% 0 and 1 are treated the same (empty)
userInputs.windows = input(['Specify number and type of spectral windows (cell array {non-negative integer, non-empty string}).',...
    newline,'(Enter {2,"logCos"} for two overlapping log cosine windows): ']);
% Check if input was non-empty:
if ~isempty(userInputs.windows)
    windowsClasses = {'cell'};
    windowsAttributes = {'vector','numel',2};
    try
        validateattributes(userInputs.windows,windowsClasses,windowsAttributes)
    catch
        disp('Error: The number and type of windows must be of the form {non-negative integer, non-empty string}. The number/type of windows has been set as empty.')
        userInputs.windows = [];    % Set erroneous input to be empty
        clear windowsClasses windowsAttributes
    end
    
    % Check if input is still non-empty (passed preceeding try):
    while ~isempty(userInputs.windows)
        
        % Check the first cell for a non-negative integer:
        windowNumberClasses = {'numeric'};
        windowsNumberAttributes = {'scalar','integer','>=',0};
        try
            validateattributes(userInputs.windows{1},windowNumberClasses,windowsNumberAttributes)
        catch
            disp('Error: The number of windows must be a non-negative integer. The number/type of windows has been set as empty.')
            userInputs.windows = [];    % Set erroneous input to be empty
            clear windowNumberClasses windowsNumberAttributes
            break
        end
        clear windowNumberClasses windowsNumberAttributes
        
        % Check if the input number of windows is 0 or 1
        if ismember(userInputs.windows{1},[0,1])
            disp('The input number of windows was either 0 or 1. The number/type of windows has been set as empty.')
            userInputs.windows = [];    % Set erroneous input to be empty
            break
        end
        
        % Check the second cell for a non-empty string:
        windowTypeClasses = {'string'};
        windowsTypeAttributes = {'nonempty','scalartext'};
        try
            validateattributes(userInputs.windows{2},windowTypeClasses,windowsTypeAttributes)
        catch
            disp('Error: The window type must be a non-empty string. The number/type of windows has been set as empty.')
            userInputs.windows = [];    % Set erroneous input to be empty
            clear windowTypeClasses windowsTypeAttributes
            break
        end
        clear windowTypeClasses windowsTypeAttributes
        
        % Check the second cell for a string with no characters:
        if strcmp(userInputs.windows{2},"")
            disp('Error: The input window type is a string without characters. The number/type of windows has been set as empty.')
            userInputs.windows = [];    % Set erroneous input to be empty
            clear windowTypeClasses windowsTypeAttributes
            break
        end
        clear windowTypeClasses windowsTypeAttributes
        break
        
    end
    clear windowsClasses windowsAttributes
end

% Set parameter selection methods:
methodsClasses = {'string'};
methodsAttributes = {'vector','nonempty'};
while isempty(userInputs.methods)
    userInputs.methods = input(['Specify parameter selection method(s) (string or array of strings).',...
        newline,'(Supported methods are ',...
        convertStringsToChars(strcat("""",join(validMethods,""", """),"""")),...
        '): ']);
    try
        validateattributes(userInputs.methods,methodsClasses,methodsAttributes)
    catch
        disp('Error: The method(s) must be input as a string or array of strings.')
        userInputs.methods = [];
    end 
end
clear methodsClasses methodsAttributes

% Set the downsampling resolutions:
userInputs.resolutions = input(['Specify downsampling resolutions (array of non-negative integers).',...
    newline,'(Supported resolutions are ',...
    convertStringsToChars(join(string(validResolutions),', ')),', or empty for full problem): ']);
if ~isempty(userInputs.resolutions)
    resolutionsClasses = {'numeric'};
    resolutionsAttributes = {'vector','integer','>=',0};
    try
        validateattributes(userInputs.resolutions,resolutionsClasses,resolutionsAttributes)
    catch
        disp('Error: The downsampling resolutions must be an array of non-negative integers. The set of resolutions has been set  as empty.')
        userInputs.resolutions = [];
    end
    
    % Remove duplicate resolutions:
    resolutions2 = unique(userInputs.resolutions);
    if numel(resolutions2) < numel(userInputs.resolutions)
        disp('Duplicate downsampling resolutions detected and removed.')
    end
    userInputs.resolutions = resolutions2;
    clear resolutionsClasses resolutionsAttributes resolutions2
end

% Shuffle images in data sets:
userInputs.shuffle = input('Do you want the images to be shuffled? (Y/N): ','s');
if ~strcmp(userInputs.shuffle,'Y') && ~strcmp(userInputs.shuffle,'N')
    disp('Error: Invalid response. The images will not be shuffled.')
    userInputs.shuffle = [];    % Empty response means no shuffling
elseif strcmp(userInputs.shuffle,'N')
    userInputs.shuffle = [];
end

% Display completion message and display contents userInputs:
disp('Constuction of the userInputs structure is complete with the following fields:')
disp(struct2table(userInputs,'AsArray',true))
