%% get_AUC.m
% This function returns the area under the curve (AUC) for the drug
% concentration vs time curve corresponding to a given dosing schedule

function auc = get_AUC(schedule)
    
    t = 0:size(schedule,2)*24;
    
    dosing_schedule = zeros(1,size(t,2));
    for j=1:size(schedule,2)
        dosing_schedule(1,((j-1)*24)+1) = dosing_schedule(1,((j-1)*24)+1) + schedule(1,j);
    end
    
%     % Full version: Layer "basis" functions, 1 for drug concentration 
%     % resulting from each day of treatment
%     funcs_vals = zeros(size(schedule,2)*24);
   
    % assuming only one day of dosing
    drug_concentration = zeros(1,size(t,2));
    received_dose = false;
    for j=1:size(dosing_schedule,2)
        if dosing_schedule(1,j) ~= 0
            received_dose = true;
            initial_dose = dosing_schedule(1,j);
            t_start = j;
        end
        
        if ~received_dose
            drug_concentration(1,j) = 0;
        else
            drug_concentration(1,j) = initial_dose.*(1/2).^(((t(1,j)+1)-t_start)/3.6124);
        end
    end
    
    % Problem: why is auc smaller for dose on day 1 than on day 3?
    
    plot(t,drug_concentration, '-o')
    axis([0 size(schedule,2)*24 0 1.25*initial_dose])
    
    auc = trapz(t,drug_concentration);
end