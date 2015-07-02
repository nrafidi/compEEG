function [dummy] = WriteSize( dataLoader, outputFile )

    Data = dataLoader();
    
    fileID = fopen(outputFile, 'w');
    for indexSize = 1:length(size(Data))
        if (indexSize > 1)
            fprintf(fileID, '\t');
        end
        fprintf(fileID, '%d', size(Data, indexSize));
    end
    
    fprintf(fileID, '\n');
    
    fclose(fileID);
    
    dummy = true;

end