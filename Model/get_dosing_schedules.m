%% get_dosing_schedules.m
function unique_perms = get_dosing_schedules(n, H, M, L)
combos = nchoosek(repelem([H M L],n),n);
permutations = [];

for i = 1:size(combos,1)
    permutations = [permutations; perms(combos(i,:))];
end

unique_perms = unique(permutations, 'rows');
end