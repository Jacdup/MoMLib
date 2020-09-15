rho = 0.4;
vertices = 40; % Number of vertices
%  Contour = [0 0 0; 1 0 0; 2 0 0];
%    Contour = [0 0 0; 0.5 0 0; 1 0 0; 1.5 0 0; 2 0 0];
%    Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 0.1; 5 0 0.2; 6 0.1 0.3; 7 0.2 0.4; 8 0.25 0.5; 9 0.2 0.6; 10 0.2 0.5; 11 0.2 0.4 ]; % Hardcode some contour
 Contour = [0 0 0; 1 0 0; 2 0 0];

%  Contour = RefineMesh(Contour,1);
[points1,quadElements1, elements1] = CylMesh(Contour,vertices,rho, cyl_def);
[revised_triangles1,triangles1] = QuadtoTri(elements1,vertices, cyl_def);
%   PlotMesh(points1,elements1,triangles1, 1);
col = [0; 1];
figure
faces = elements1(1:end-1,1:2)';
faces = [faces;elements1(1:end-1,3:4)'];
patch('Faces',faces,'Vertices',points1,'FaceVertexCData',col,'FaceColor','white')
%    trimesh(quadElements1,points1(:,1),points1(:,2),points1(:,3)) ;
%   [node_coords,quad_Elements] = MeshRectanglularPlate(1,1,4,4) ;

%   Delta = importDelta('C:\Users\19083688\Desktop\Delta.txt', 1, 1500);
%Delta = importDelta('C:\Users\19083688\Desktop\MoM Codes\V_3.6_feat_plane_wave_MBF_qmom\Delta.txt', 1, 1500);
% ruv = importruv('C:\Users\19083688\Desktop\ruv.txt', 1, 1500);
 Delta = importDeltaNew('C:\Users\19083688\Desktop\Masters\MoMLib\Delta.txt',[1, Inf]);
 Ruv = importruv("C:\Users\19083688\Desktop\Masters\MoMLib\Ruv.txt", [1, Inf]);
 Delta1 = importDeltaNew('C:\Users\19083688\Desktop\Masters\MoMLib\Delta1.txt',[1, Inf]);
 Ruv1 = importruv("C:\Users\19083688\Desktop\Masters\MoMLib\Ruv1.txt", [1, Inf]);
%   quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), Delta(:,1)+Delta1(:,1), Delta(:,2)+Delta1(:,2), Delta(:,3)+Delta1(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
%    quiver3(Ruv1(1:26,1),Ruv1(1:26,2),Ruv1(1:26,3), Delta(1:26,1), Delta(1:26,1), Delta(1:26,1), 1,'LineWidth', 2,'MaxHeadSize',5);
%  Delta2 = Delta1 + Delta;
%figure;
%h = subplot(1,3,1);
  hold on
  for i = 1:4
  domain = i:64:length(Ruv)/2;
  domain2 = length(Ruv)/2 + i:32:length(Ruv);
%     domain3 = 3:16:length(Ruv);
% quiver3(ruv_large(:,1),ruv_large(:,2),ruv_large(:,3), DeltaLarge(:,1), DeltaLarge(:,2), DeltaLarge(:,3), 5);
quiver3(Ruv(domain,1),Ruv(domain,2),Ruv(domain,3), Delta(domain,1), Delta(domain,2), Delta(domain,3), 1.5,'LineWidth', 1.5,'MaxHeadSize',3);
quiver3(Ruv(domain2,1),Ruv(domain2,2),Ruv(domain2,3), Delta(domain2,1), Delta(domain2,2), Delta(domain2,3), 1.5,'LineWidth', 1.5,'MaxHeadSize',3);
  end 
  axis equal
% % Delta(:,1) = 0;
% % Delta1(:,1) = 0;
% % First edge (DOFS 1+2)
% rho = [-max(Delta(1:26,2)),-max(Delta(1:26,2));-max(Delta1(1:26,2)),max(Delta1(1:26,2))];
% B1   = [1;0];
% B1 = B1.* [sind(45)';sind(45)'];
% X1   = rho\B1;
% 
% % Straight edge (DOFS 3+4)
% rho = [-max(Delta(27:52,2)),-max(Delta(27:52,2));max(Delta1(27:52,2)),-max(Delta1(27:52,2))];
% B2   = [1;1];
% X2   = rho\B2;
% 
% % Second edge (DOFS 5+6)
% rho = [-max(Delta(53:end,2)),-max(Delta(53:end,2));max(Delta1(53:end,2)),-max(Delta1(53:end,2))];
% B3   = [0;1];
% B3 = B3.* [sind(45)';sind(45)'];
% X3   = rho\B3;
% 
% RWG_Dofs_1 = Delta1(1:26,:);
% RWG_Dofs_2 = Delta1(27:52,:);
% RWG_Dofs_3 = Delta1(53:end,:);
% Lin_Dofs_1 = Delta(1:26,:);
% Lin_Dofs_2 = Delta(27:52,:);
% Lin_Dofs_3 = Delta(53:end,:);
% 
% R1 = (X1(1) * RWG_Dofs_1) + (X1(2) * Lin_Dofs_1);
% R2 = (X2(1) * RWG_Dofs_2) + (X2(2) * Lin_Dofs_2);
% R3 = (X3(1) * RWG_Dofs_3) + (X3(2) * Lin_Dofs_3);
% % R1(:,1) = 0;
% % R2(:,1) = 0;
% % R3(:,1) = 0;
% 
% % R = R1 + R2 + R3;
% 
% %  [node_coords,quad_Elements] = MeshRectanglularPlate(a(iter),a(iter),8,8) ;
% %  hold on
% % quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), X(1)*Delta(:,1), X(1)*Delta(:,2), X(1)*Delta(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), X(2)*Delta1(:,1), X(2)*Delta1(:,2), X(2)*Delta1(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% %  [node_coords,quad_Elements] = MeshRectanglularPlate(a(iter),a(iter),8,8) ;
% %  hold on
% % %  quiver3(Ruv(1:26,1),Ruv(1:26,2),Ruv(1:26,3), R(:,1), R(:,2), R(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% %  quiver3(Ruv(1:26,1),Ruv(1:26,2),Ruv(1:26,3), R1(:,1), R1(:,2), R1(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% %   quiver3(Ruv(27:52,1),Ruv(27:52,2),Ruv(27:52,3), R2(:,1), R2(:,2), R2(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% %     quiver3(Ruv(53:end,1),Ruv(53:end,2),Ruv(53:end,3), R3(:,1), R3(:,2), R3(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % % quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), Delta2(:,1), Delta2(:,2), Delta2(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % % quiver3(Edge3ruv(:,1),Edge3ruv(:,2),Edge3ruv(:,3), Edge1(:,1), Edge1(:,2), Edge1(:,3), 0.5);
% % axis equal
% % triplot(elem_nodes,node_coords(:,1),node_coords(:,2),'red','LineWidth',linewidth_triplot);
% % hold on;
% % quiver(coord_xy(:,1),coord_xy(:,2),bf1(:,1),bf1(:,2),'LineWidth',linewidth_quiver);
% % axis(axis_limits);
% % set(gca,'DataAspectRatio',[1 1 1]); % axes equal;