function table = pdReadRawCSV(filename,n)
    if n == 1; varName = "IN";
    elseif n == 2; varName = "B";
    elseif n == 3; varName = "D";
    elseif n == 4; varName = "NP";
    end
    %--- Setup the Import Options and import the Data
    opts = delimitedTextImportOptions("NumVariables",4);
    % Specify range and delimiter
    opts.DataLines = [2,Inf];
    opts.Delimiter = ",";
    % Specify column names and types
    opts.VariableNames = ["IN","B","D","NP"]; % IN = "ki,koff,kon,L"
    opts.SelectedVariableNames = varName;
    opts.VariableTypes = ["string","string","string","string"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts,["IN","B","D","NP"],"WhitespaceRule","preserve");
    opts = setvaropts(opts,["IN","B","D","NP"],"EmptyFieldRule","auto");
    %--- Import the Data
    table = readtable(filename,opts);
end