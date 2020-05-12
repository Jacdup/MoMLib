function [x_axis, current_norm, current_MBF] = plotNodeCurrent(I_vec,I_vec_norm,DOF_mat1,DOF_mat2,DOF_mat3,DOF_mat, numNodes, endcap)
iter = 0;
% Current magnitude at all nodes
for node = 1:numNodes
    iter = iter + 1;
    node_to_plot = 4;
    dofs = 1:2:length(DOF_mat1);
    dofsLin = 2:2:length(DOF_mat1);
    
    dofs_circ = 1:2:length(DOF_mat);
    dofs_circ_lin = 2:2:length(DOF_mat);
%     x_axis = 0:360/vertices:359;
    x_axis = linspace(0,360,length(dofs));
    x_axis_circ = linspace(0,360,length(dofs_circ));
    if (endcap) % Exclude zeros
        DOF_mat2 = DOF_mat2(:,2:end);
        DOF_mat3 = DOF_mat3(:,1:end-1);
    end
%     current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node)) + I_vec_norm(DOF_mat1(dofsLin,node)));
%      current_MBF(:,iter) = abs(I_vec(DOF_mat1(dofs,node)) + I_vec(DOF_mat1(dofsLin,node)));
    current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node)) + I_vec_norm(DOF_mat1(dofsLin,node)) + I_vec_norm(DOF_mat2(dofs,node)) + I_vec_norm(DOF_mat2(dofsLin,node)) + I_vec_norm(DOF_mat3(dofs,node)) + I_vec_norm(DOF_mat3(dofsLin,node))) ;
    current_MBF(:,iter) = abs(I_vec(DOF_mat1(dofs,node)) + I_vec(DOF_mat1(dofsLin,node)) + I_vec(DOF_mat2(dofs,node)) + I_vec(DOF_mat2(dofsLin,node)) + I_vec(DOF_mat3(dofs,node)) + I_vec(DOF_mat3(dofsLin,node))) ;

%     current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node))+I_vec_norm(DOF_mat2(dofs,node))+I_vec_norm(DOF_mat3(dofs,node)));
%     current_MBF(:,iter)  = abs(I_vec(DOF_mat1(dofs,node))+I_vec(DOF_mat2(dofs,node))+I_vec(DOF_mat3(dofs,node)));
    
%      current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node)));
%     current_MBF(:,iter)  = abs(I_vec(DOF_mat1(dofs,node)));
    
    current_norm_circ(:,iter) = abs(I_vec_norm(DOF_mat(dofs_circ,node)) + I_vec_norm(DOF_mat(dofs_circ_lin,node)));
    current_MBF_circ(:,iter) = abs(I_vec(DOF_mat(dofs_circ,node))+ I_vec(DOF_mat(dofs_circ_lin,node)));

%     current_norm_circ(:,iter) = abs(I_vec_norm(DOF_mat(dofs_circ,node)));
%     current_MBF_circ(:,iter) = abs(I_vec(DOF_mat(dofs_circ,node)));
    
end
    iter = iter + 1;
%      current_norm_circ(:,iter) = abs(I_vec_norm(DOF_mat(dofs_circ,node+1)));
%     current_MBF_circ(:,iter) = abs(I_vec(DOF_mat(dofs_circ,node+1)));

     current_norm_circ(:,iter) = abs(I_vec_norm(DOF_mat(dofs_circ,node+1)) + I_vec_norm(DOF_mat(dofs_circ_lin,node+1)));
    current_MBF_circ(:,iter) = abs(I_vec(DOF_mat(dofs_circ,node+1)) + I_vec(DOF_mat(dofs_circ_lin,node+1)));


    figure
%         plot(x_axis_circ,current_norm_circ);
%         hold on
%         plot(x_axis_circ,current_MBF_circ);
%         hold off
surf(x_axis, linspace(0,2,numNodes), (current_norm - current_MBF)');
xlabel('Circumference (\phi)');
ylabel('Cylinder length (m)');
zlabel('Axial Current magnitude');

figure
surf(x_axis_circ, linspace(0,2,numNodes+1), (current_norm_circ - current_MBF_circ )');
xlabel('Circumference (\phi)');
ylabel('Cylinder length (m)');
zlabel('Azimuthal Current magnitude');
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
    