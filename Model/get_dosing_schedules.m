%% get_dosing_schedules.m
% DESCRIPTION:
% This function creates and returns a matrix of all the possible permutations
% of dosing schedules with three doses (high=H, medium=M, low=L) over n days. 

% AUTHORS: Zach Schlamowitz and Andrew Paek, 10/19/21

function treatments = get_dosing_schedules(n, H, M, L)
        treatments = zeros(3^n, n);
        
        for j = 1:n  % cols
            for i = 1:3^n  % rows
                if  mod(fix((i-1)/(3^(j-1))),3)    == 0
                    treatments(i,n+1-j) = H;
                elseif mod(fix((i-1)/(3^(j-1))),3) == 1
                    treatments(i,n+1-j) = M;
                elseif mod(fix((i-1)/(3^(j-1))),3) == 2
                    treatments(i,n+1-j) = L;
                end
                
                % disp(treatments)
            end
        end
end

% OLD VERSION:
% -------------
% function unique_perms = get_dosing_schedules(n, H, M, L)
% % Identify all possible combinations of n choices of H/M/L; initialize
% % matrix for permutations of these
% combos = nchoosek(repelem([H M L],n),n);
% permutations = [];
% 
% % Get all permutations of each combination
% for i = 1:size(combos,1)
%     permutations = [permutations; perms(combos(i,:))];
% end
% 
% % Remove duplicates
% unique_perms = unique(permutations, 'rows');
% 
% end

% get_dosing_scheduless(7, 3,2,1)
