function dataname = loadData()
% This function prompts the user to select which data file to load. If no
% such data file exists, a message is printed to the user indicating that 
% no data name is returned. If the specified data file exists, the output
% dataname is the name of the data file. Use the 'load' command to load the
% data contained in the file named dataname.

% Ask the user if they want to load data:
specs = questdlg('Would you like to load data?',...
    'Load data?','Yes','No','No');

switch specs
    case 'No'
        dataname = [];
        disp('No data has been loaded.')
    case 'Yes'
        prompt = {'Enter dimension (1D, 2D, or 3D):',...
            'Enter test function number (1, 2, or 3):',...
            'Enter SNR:',...
            'Enter width parameter of Gaussian kernel:',...
            'Enter number of noise realizations:'};
        title = 'Data specifications';
        specs = inputdlg(prompt,title);
        specs{3} = num2str(str2double(specs{3}),'%02.f');
        specs{4} = num2str(str2double(specs{4}),'%03.f');
        specs{5} = num2str(str2double(specs{5}),'%02.f');
        dataname = ['Data',specs{1},'_F',specs{2},'_S',specs{3},'_W',...
            specs{4},'_R',specs{5},'.mat'];
        
        % Check to see if the file exists
        if exist(dataname,'file') == 0  
            disp(['Error, no such data file exists yet: ' dataname])
            dataname = 'None';  % Delete invalid data name
        end
        
end

end
