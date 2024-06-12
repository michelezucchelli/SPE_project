function config = load_config(file_path)    
    % sets options
    opts = delimitedTextImportOptions('CommentStyle','#','delimiter',...
        ': ','VariableNamingRule','preserve');
    
    % reads in config as a table, rotates its orientation, converts it
    % to a struct, and removes extra row storing original variable names
    config_struct = rmfield(table2struct(rows2vars(readtable(file_path,...
        opts))),'OriginalVariableNames');
    
    % field names stored in array
    fields = fieldnames(config_struct);
    
    % structs storing key-value pairs
    keys = config_struct(1);
    values = config_struct(2);
    
    % initializes config dictionary to return
    config = dictionary;
    
    % simplifies struct produces by table2struct, also converting all
    % fields from chars to strings
    for i = 1:length(fields)        
        % extracts key and value
        key = keys.(fields{i});
        value = values.(fields{i});        
        % converts any string values to char arrays
        value = convertStringsToChars(value); 


        % empty config
        if isempty(value)
            if key(end) == ':'
                key = key(1:end-1);
            end
            value = [];            
        % % logical config
        % elseif strcmpi(value,'false')
        %     value = false;        
        % elseif strcmpi(value,'true')
        %     value = true;            
        % array config
        elif strcmpi(value(1),'[')           
            % removes any extra whitespace
            value = strrep(value,', ',',');            
            % removes brackets and ensures the array is a string array
            value = string(split(value(2:end-1),',')).';            
            % converts to double or logical array if value elements are 
            % numeric or logical
            if isnumeric(value(1)) || islogical(value(1))
                value = str2num(convertStringsToChars(strjoin(value)));
            end            
        end

        config{key} = value;        
    end
    
end
%% Function Description
% Loads a simple configuration file into a dictionary.
% Input: 
%    - filepath: path of the file ...
% Output: 
%    - config: dictionary in the form key-values
% 
% This function is based on the one written by Tamas Kis 
% https://it.mathworks.com/matlabcentral/fileexchange/129884-simple-
% configuration-file-format-for-matlab
