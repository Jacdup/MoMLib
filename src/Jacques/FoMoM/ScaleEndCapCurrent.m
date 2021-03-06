function [I_vec] = ScaleEndCapCurrent(DOF_mat1, DOF_mat,theta1, theta2, I_vec, cyl_def)

DOF_mat(:,~any(DOF_mat,1)) = []; % Remove zero columns
DOF_mat1(:,~any(DOF_mat1,1)) = []; % Remove zero columns
if (cyl_def.firstNode == "endCap")
    first_endcap_axial       = DOF_mat1(:,1);
    if cyl_def.lastNode == "endCap"
        first_endcap_circ        = nonzeros(DOF_mat(:,end-1));
    else
        first_endcap_circ        = nonzeros(DOF_mat(:,end));
    end
    currents_endcap1         = I_vec(first_endcap_axial);
    current_maxVal_endcap1   = max(currents_endcap1);
    I_vec(first_endcap_circ(2:2:end)) = dot(I_vec(first_endcap_circ(2:2:end)),repelem((current_maxVal_endcap1), length(first_endcap_circ)/2)); % So that all values are smaller than maxVal
end

if cyl_def.lastNode == "endCap"
    second_endcap_axial      = DOF_mat1(:,end);
    second_endcap_circ       = nonzeros(DOF_mat(:,end)); % Remember to make this just 'end' when not conn
    currents_endcap2         = I_vec(second_endcap_axial);
    current_maxVal_endcap2   = max(currents_endcap2);
    I_vec(second_endcap_circ(2:2:end)) = dot(I_vec(second_endcap_circ(2:2:end)),repelem((current_maxVal_endcap2), length(second_endcap_circ)/2));
end

% I_vec(first_endcap_circ) = I_vec(first_endcap_circ)*real(current_maxVal_endcap1); % So that all values are smaller than maxVal
% I_vec(second_endcap_circ) = I_vec(second_endcap_circ)*real(current_maxVal_endcap2);
% I_vec(first_endcap_circ(1:2:end)) = I_vec(first_endcap_circ(1:2:end)).*-1; % So that all values are smaller than maxVal
% I_vec(second_endcap_circ(1:2:end)) = I_vec(second_endcap_circ(1:2:end)).*-1;

% I_vec(first_endcap_circ(2:2:end)) = I_vec(first_endcap_circ(2:2:end))*real(current_maxVal_endcap1); % So that all values are smaller than maxVal
% I_vec(second_endcap_circ(2:2:end)) = I_vec(second_endcap_circ(2:2:end))*real(current_maxVal_endcap2);



end