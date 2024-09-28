function [outputArg1,outputArg2] = table_header(filename, method_names)
fileID = fopen(filename, 'w'); 
n_methods = length(method_names);
fprintf(fileID, '\\begin{tabular}{|l|l|');
for i =1:n_methods
    fprintf(fileID, 'c');
end
fprintf(fileID, '|}\n \\hline \n Scene & Metric & ');
for method_id = 1:n_methods
    fprintf(fileID, strrep(method_names(method_id),'_', ' '));
    if method_id < n_methods
        fprintf(fileID, ' & ');
    end
end
fprintf(fileID, '\\\\ \n \\hline \n');
end

