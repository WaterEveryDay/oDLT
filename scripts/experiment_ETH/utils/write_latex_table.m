function write_latex_table(scene, method_names, ...
    aggregate_err_rot_all, aggregate_err_pos_all ...
    , aggregate_err_reproj_all, aggregate_time_all, n_avg, filename)
% Open the file in append mode
fileID = fopen(filename, 'a'); % 'a' mode for appending

% Number of methods
n_methods = length(method_names);

% Function to format a row with bold and underline
    function write_row(label, values, n_avg, row_number)
% Find the indices of the lowest and second lowest scores
[sorted_values, sorted_indices] = sort(values);
lowest_idx = sorted_indices(1);
if n_methods > 1
    second_lowest_idx = sorted_indices(2);
    third_lowest_idx = sorted_indices(3);
else
    second_lowest_idx = lowest_idx; % Handle case with only one method
    third_lowest_idx = lowest_idx; % Handle case with only one method
end

% if first row, then print scene name
if row_number == 1
    fprintf(fileID, '%s & %s & ', strrep(scene, '_', ' ') , label);
elseif row_number == 2
    fprintf(fileID, 'average $n$ %d & %s & ', floor(n_avg) , label); 
else
    fprintf(fileID, ' & %s & ', label); 
end

for method_id = 1:n_methods
    if method_id == lowest_idx
        fprintf(fileID, '\\textbf{%8.5f}', values(method_id));
    elseif method_id == second_lowest_idx
        fprintf(fileID, '\\underline{%8.5f}', values(method_id));
    elseif method_id == third_lowest_idx
        fprintf(fileID, '\\textit{%8.5f}', values(method_id));
    else
        fprintf(fileID, '%8.5f', values(method_id));
    end
    if method_id < n_methods
        fprintf(fileID, ' & ');
    end
end
fprintf(fileID, ' \\\\ \n');
end


% Write Rotation RMSE section
write_row('Rot. RMSE (deg)', aggregate_err_rot_all(1, :), n_avg, 1);

% Write Position RMSE section
write_row('Pos. RMSE (m)', aggregate_err_pos_all(1, :), n_avg, 2);

% Write Mean Geometric Error section
write_row('Mean Reproj. Err. (pixel)', aggregate_err_reproj_all(1, :), n_avg, 3);

% Write Mean Run Time section
write_row('Mean Run Time (ms)', 1000 * aggregate_time_all(1, :), n_avg, 4);

% End the table with a horizontal line
fprintf(fileID, '\\hline\n');

% Close the file
fclose(fileID);
end

