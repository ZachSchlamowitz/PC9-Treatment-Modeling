% For 1uM Gefitnib:
num_divs_per_cell = sum(Clean_Gef_1uM_200cells.divisionMatrixDataset,2);
div_flags = zeros(numel(num_divs_per_cell),1);
for i=1:numel(num_divs_per_cell)
    if mod(num_divs_per_cell(i),2) ~= 0
        fprintf('flag cell ')
        fprintf(num2str(i))
        fprintf('\n')
        div_flags(i) = 1;
    end
end

total_divs = sum(num_divs_per_cell);
for i=[1 18 22 26 41 48 66 125 126]
    total_divs = total_divs - num_divs_per_cell(i);
end

gef_stats_1uM.div_events_per_cell = num_divs_per_cell;
gef_stats_1uM.div_flags = div_flags;
gef_stats_1uM.num_divs = total_divs/2;
gef_stats_1uM.prob_div = gef_stats_1uM.num_divs/gef_stats_1uM.num_cells_good4divs

% For 0.1uM Gefitnib:
num_divs_per_cell = sum(CleanGef0_1.divisionMatrixDataset,2);
div_flags = zeros(numel(num_divs_per_cell),1);
for i=1:numel(num_divs_per_cell)
    if mod(num_divs_per_cell(i),2) ~= 0
        fprintf('flag cell ')
        fprintf(num2str(i))
        fprintf('\n')
        div_flags(i) = 1;
    end
end

total_divs = sum(num_divs_per_cell);
for i=[33 67 68 134 165]
    total_divs = total_divs - num_divs_per_cell(i);
end

gef_stats_0p1uM.div_events_per_cell = num_divs_per_cell;
gef_stats_0p1uM.div_flags = div_flags;
gef_stats_0p1uM.num_divs = total_divs/2;
gef_stats_0p1uM.prob_div = gef_stats_0p1uM.num_divs/gef_stats_0p1uM.num_cells_good4divs
