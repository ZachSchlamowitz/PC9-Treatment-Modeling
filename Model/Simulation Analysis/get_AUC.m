%% get_AUC.m
% This function returns the area under the curve (AUC) for the drug
% concentration vs time curve corresponding to a given dosing schedule

function auc = get_AUC(schedule)
    
    drug_half_life = 3.6124; % (hrs) default: 48; 99% decay in 24 hrs: 3.6124

    t = 0:size(schedule,2)*24;
    
    dosing_schedule = zeros(1,size(t,2));
    for j=1:size(schedule,2)
        dosing_schedule(1,((j-1)*24)+1) = dosing_schedule(1,((j-1)*24)+1) + schedule(1,j);
    end
    
%     % Full version: Layer "basis" functions, 1 for drug concentration 
%     % resulting from each day of treatment
%     funcs_vals = zeros(size(schedule,2)*24);


    drug_concentrations = zeros(size(schedule,2),size(t,2));
    for i = 1:size(schedule,2)            
        initial_dose = schedule(1,i);
        t_start = 1+(i-1)*24;
        
        for j=1:size(dosing_schedule,2)
            if j < t_start
                drug_concentrations(i,j) = 0;
            else
                drug_concentrations(i,j) = initial_dose.*(1/2).^(((t(1,j)+1)-t_start)/drug_half_life);
            end
        end
        
        
%         % assuming only one day of dosing
%         received_dose = false;
%         for j=1:size(dosing_schedule,2)
%             if dosing_schedule(1,j) ~= 0
%                 received_dose = true;
%                 initial_dose = dosing_schedule(1,j);
%                 t_start = j;
%             end
% 
%             if ~received_dose
%                 drug_concentration(i,j) = 0;
%             else
%                 drug_concentration(i,j) = initial_dose.*(1/2).^(((t(1,j)+1)-t_start)/48);  % 99% in 24hrs half life: 3.6124
%             end
%         end
    end
    
    drug_concentration = sum(drug_concentrations,1);
    % Problem: why is auc smaller for dose on day 1 than on day 3?
    % Answer: It's the way that the curve is drawn. When you start on day
    % 1, the curve starts high, so the left boundary of area is a vertical
    % line. However, when you start on a day in the middle, you go from
    % say, (47,0) to (48,10), which creates an additional triangle of area
    % (1/2)*dose
    
    plot(t,drug_concentration, '-o')
    axis([0 150 0 1.2])
    
    auc = trapz(t,drug_concentration);
end