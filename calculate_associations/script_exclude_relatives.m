    %% Exclude relatives from the analysis
    % The script reads a list of first-degree relatives from a file and
    % compares it with the subject IDs. If a subject ID is found in the
    % list, the corresponding columns in the data matrices are set to NaN
    % (and excluded from the analyses).

    % Read the list of first-degree relatives from a file
    firstDegreeList = csvread(fullfile(absolute_path, 'data', ...
        'relatedSubjects.txt'));

    % Convert subject IDs to numbers
    subjectID_num = cellfun(@str2num, subjectID);

    % Find the indices that are in the first-degree relatives list
    indx = ismember(subjectID_num, firstDegreeList);    

    % Set the corresponding columns in the data matrices to NaN
    surfaces(:, indx) = NaN;
    thicknesses(:, indx) = NaN;
    volumes(:, indx) = NaN;
    connvF(:, indx) = NaN;
    connvS(:, indx) = NaN;
    tfmri(:, indx) = NaN;

    % Present the numbers
    fprintf('%g relatives excluded from %g total participants.\n', ...
        nnz(indx), length(tfmri));