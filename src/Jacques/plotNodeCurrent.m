function [x_axis, current_norm, current_MBF] = plotNodeCurrent(I_vec,I_vec_norm,DOF_mat1,DOF_mat2,DOF_mat3,DOF_mat, numNodes)
iter = 0;
% Current magnitude at all nodes
for node = 1:numNodes
    iter = iter + 1;
    node_to_plot = 4;
    dofs = 1:2:length(DOF_mat1);
    dofs_circ = 1:4:length(DOF_mat);
%     x_axis = 0:360/vertices:359;
    x_axis = linspace(0,360,length(dofs));
    x_axis_circ = linspace(0,360,length(dofs_circ));

    current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node))+I_vec_norm(DOF_mat2(dofs,node))+I_vec_norm(DOF_mat3(dofs,node)));
    current_MBF(:,iter)  = abs(I_vec(DOF_mat1(dofs,node))+I_vec(DOF_mat2(dofs,node))+I_vec(DOF_mat3(dofs,node)));
    
    current_norm_circ(:,iter) = abs(I_vec_norm(DOF_mat(dofs_circ,node)));
    current_MBF_circ(:,iter) = abs(I_vec(DOF_mat(dofs_circ,node)));
    

        
end

    figure
%         plot(x_axis_circ,current_norm_circ);
%         hold on
%         plot(x_axis_circ,current_MBF_circ);
%         hold off
surf(x_axis_circ, linspace(1,2,numNodes), (current_norm_circ)');
% hold on
% surf(x_axis, 1:numNodes, current_MBF');
% hold off
% axis equal
%         plot(x_axis,current_norm);
%         hold on
%         plot(x_axis,current_MBF);
%         hold off
%         legend('Real current magnitude', 'MBF current magnitude');
end
    