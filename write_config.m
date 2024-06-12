function [] = write_config(filename, key, values)
    % if the file already exists, then open it in
    if isfile(filename)
        % open the file
        fileId = fopen(filename, 'a');
        % load the file
        fileL = load_config(filename);
        % if the key is already present in the file loaded, then do not
        % change the associated values 
        try
           val = fileL(num2str(key));
        catch % otherwise add the new key-value pair
           % disp('Key not already present, then append it to the file');
            fprintf(fileId, '%d: ', key);
            fprintf(fileId, '%f, ', values);
            fprintf(fileId, '\n');
        end 
    else % otherwise create the file and write the key-value pair
        fileId = fopen(filename, 'w');
        fprintf(fileId, '%d: ', key);
        fprintf(fileId, '%f, ', values);
        fprintf(fileId, '\n');
    end
    fclose(fileId);
end

%% Function Description
% Write on a file a key-values pair (dictionary)
% Inputs:
%    - filename: name of the file ...
%    - key: unique
%    - values: array of values associated to a key