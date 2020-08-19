function [x_axis, current_norm, current_MBF] = plotNodeCurrent(I_vec,I_vec_norm,DOF_mat1,DOF_mat2,DOF_mat3,DOF_mat, numNodes, cyl_def)
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
    if (cyl_def.firstNode == "endcap" && cyl_def.lastNode == "endcap" ) % Exclude zeros
        DOF_mat2 = DOF_mat2(:,2:end);
        DOF_mat3 = DOF_mat3(:,1:end-1);
    end
   
    %     current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node)) + I_vec_norm(DOF_mat1(dofsLin,node)));
%      current_MBF(:,iter) = abs(I_vec(DOF_mat1(dofs,node)) + I_vec(DOF_mat1(dofsLin,node)));
  current_norm(:,iter) = real(I_vec_norm(DOF_mat1(dofs,node)) + I_vec_norm(DOF_mat2(dofs,node)) + I_vec_norm(DOF_mat3(dofs,node)) ) ;
%     current_norm(:,iter) = abs(I_vec_norm(DOF_mat1(dofs,node)) + I_vec_norm(DOF_mat1(dofsLin,node)) + I_vec_norm(DOF_mat2(dofs,node)) + I_vec_norm(DOF_mat2(dofsLin,node)) + I_vec_norm(DOF_mat3(dofs,node)) + I_vec_norm(DOF_mat3(dofsLin,node))) ;
    current_MBF(:,iter) = real(I_vec(DOF_mat1(dofs,node)) + I_vec(DOF_mat1(dofsLin,node)) + I_vec(DOF_mat2(dofs,node)) + I_vec(DOF_mat2(dofsLin,node)) + I_vec(DOF_mat3(dofs,node)) + I_vec(DOF_mat3(dofsLin,node))) ;

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

% figure('DefaultAxesFontSize',14);
% surf(x_axis, linspace(0,2,numNodes), (current_MBF)');
% xlabel('$\phi$ (deg)');
% ylabel('($z/\lambda$)');
% zlabel('$\textrm{Re}(J_z^{MBF})$');
% % zlabel('$|J_z^{MBF}|$');
% 
% figure('DefaultAxesFontSize',14);
% surf(x_axis, linspace(0,2,numNodes), (current_norm)');
% xlabel('$\phi$ (deg)');
% ylabel('($z/\lambda$)');
% zlabel('$\textrm{Re}(J_z^{RWG})$');
% % zlabel('$|J_z^{RWG}|$');
% 
% figure('DefaultAxesFontSize',14);
% surf(x_axis_circ, linspace(0,2,numNodes+1), (current_MBF_circ )');
% xlabel('$\phi$ (deg)');
% ylabel('Cylinder length ($z/\lambda$)');
% zlabel('$|J_\phi^{MBF}|$');

for node_to_plot = 11:11
     current_MBF_node(:) = imag((I_vec(DOF_mat1(dofs,node_to_plot)) + I_vec(DOF_mat1(dofsLin,node_to_plot)) + I_vec(DOF_mat2(dofs,node_to_plot)) + I_vec(DOF_mat2(dofsLin,node_to_plot)) + I_vec(DOF_mat3(dofs,node_to_plot)) + I_vec(DOF_mat3(dofsLin,node_to_plot)))) ;
    current_norm_node(:) = imag(I_vec_norm(DOF_mat1(dofs,node_to_plot)) + I_vec_norm(DOF_mat2(dofs,node_to_plot)) + I_vec_norm(DOF_mat3(dofs,node_to_plot)) ) ;
    
figure('DefaultAxesFontSize',14);
plot(x_axis,current_MBF_node);
hold on
plot(x_axis,current_norm_node);

ylabel('$$\textrm{Re}(J_z)$$');
xlabel('$\phi$ (deg)');
% ylabel('$|J_z|$');
legend('MBF', 'RWG');
end
end
    