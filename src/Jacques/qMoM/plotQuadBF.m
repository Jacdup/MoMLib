rho = 0.2;
vertices = 8; % Number of vertices
%  Contour = [0 0 0; 1 0 0; 2 0 0];
%    Contour = [0 0 0; 0.5 0 0; 1 0 0; 1.5 0 0; 2 0 0];
%    Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 0.1; 5 0 0.2; 6 0.1 0.3; 7 0.2 0.4; 8 0.25 0.5; 9 0.2 0.6; 10 0.2 0.5; 11 0.2 0.4 ]; % Hardcode some contour
 Contour = [0 0 0; 0.5 0 0; 1 0 0];

%  Contour = RefineMesh(Contour,1);
[points1,quadElements1, elements1] = CylMesh(Contour,vertices,rho, cyl_def);
[revised_triangles1,triangles1] = QuadtoTri(elements1,vertices, cyl_def);
%   PlotMesh(points1,elements1,triangles1, 1);
col = [0; 1];
% figure
% faces = elements1(1:end-1,1:2)';
% faces = [faces;elements1(1:end-1,3:4)'];
% patch('Faces',faces,'Vertices',points1,'FaceVertexCData',col,'FaceColor','white')
%    trimesh(quadElements1,points1(:,1),points1(:,2),points1(:,3)) ;
  [node_coords,quad_Elements] = MeshRectanglularPlate(1,1,2,2) ;

%   Delta = importDelta('C:\Users\19083688\Desktop\Delta.txt', 1, 1500);
%Delta = importDelta('C:\Users\19083688\Desktop\MoM Codes\V_3.6_feat_plane_wave_MBF_qmom\Delta.txt', 1, 1500);
% ruv = importruv('C:\Users\19083688\Desktop\ruv.txt', 1, 1500);
 Delta = importDeltaNew('C:\Users\19083688\Desktop\Masters\MoMLib\Delta.txt',[1, Inf]);
 Ruv = importruv("C:\Users\19083688\Desktop\Masters\MoMLib\Ruv.txt", [1, Inf]);
 Delta1 = importDeltaNew('C:\Users\19083688\Desktop\Masters\MoMLib\Delta1.txt',[1, Inf]);
 Ruv1 = importruv("C:\Users\19083688\Desktop\Masters\MoMLib\Ruv1.txt", [1, Inf]);
   hold on 
%    Delta(:,:) = Delta(:,:) * 
%  quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), Delta(:,1)+Delta1(:,1), Delta(:,2)+Delta1(:,2), Delta(:,3)+Delta1(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
%  r1 = 1:25
% Ruv1(:,3) = 0;
% Ruv1(:,1) = Ruv1(:,1)*0.5;
% Ruv1(26:50,:) = [];
% Ruv1(101:125,:) = [];
 Ruv1 = sortrows(Ruv1,2);
% Ruv1 = unique(Ruv1,'rows');
% Ruv1(24:2:76,:) = [];
% figure
% Delta(:,[1 2]) = Delta(:,[2 1]);
%  Delta(1:length(Delta)/2,:) = sortrows(Delta(1:length(Delta)/2,:), 'descend');
%   Delta(length(Delta)/2:length(Delta),:) = sortrows(Delta(length(Delta)/2:length(Delta),:));
quiver3(Ruv1(1:2:end,1),Ruv1(1:2:end,2),Ruv1(1:2:end,3), Delta(:,1), Delta(:,2), Delta(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
hold on
quiver3(Ruv1(2:2:end,1),Ruv1(2:2:end,2),Ruv1(2:2:end,3), Delta(:,1), Delta(:,2), Delta(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
axis equal
%  Delta2 = Delta1 + Delta;
%figure;
%h = subplot(1,3,1);

%   for i = 1:4
%   domain = i:64:length(Ruv)/2;
%   domain2 = length(Ruv)/2 + i:32:length(Ruv);
% %     domain3 = 3:16:length(Ruv);
% % quiver3(ruv_large(:,1),ruv_large(:,2),ruv_large(:,3), DeltaLarge(:,1), DeltaLarge(:,2), DeltaLarge(:,3), 5);
% quiver3(Ruv(domain,1),Ruv(domain,2),Ruv(domain,3), Delta(domain,1), Delta(domain,2), Delta(domain,3), 1.5,'LineWidth', 1.5,'MaxHeadSize',3);
% quiver3(Ruv(domain2,1),Ruv(domain2,2),Ruv(domain2,3), Delta(domain2,1), Delta(domain2,2), Delta(domain2,3), 1.5,'LineWidth', 1.5,'MaxHeadSize',3);
%   end 
%   axis equal
% Delta(:,1) = 0;
% Delta1(:,1) = 0;
lim1 = 1:50;
lim1_1 = 26:50;
lim2 = 50:100;
lim2_2 = 100:125;
% Ruv2_lim = 51:100;
% Ruv2_lim = 26:2:125;
lim3 = 101:150;
% First edge (DOFS 1+2)
rho = [-max(Delta(lim1,2)),-max(Delta(lim1,2));-max(Delta1(lim1,2)),max(Delta1(lim1,2))];
B1   = [1;0];
B1 = B1.* [sind(45)';sind(45)'];
X1   = rho\B1;

% Straight edge (DOFS 3+4)
rho = [1,1;1,-1];
B2   = [1;1];
X2   = rho\B2;

% Second edge (DOFS 5+6)
rho = [-max(Delta(lim3,2)),-max(Delta(lim3,2));max(Delta1(lim3,2)),-max(Delta1(lim3,2))];
B3   = [0;-1];
B3 = B3.* [sind(45)';sind(45)'];
X3   = rho\B3;

RWG_Dofs_1 = Delta1(lim1,:);
RWG_Dofs_2 = Delta1(lim2,:);
RWG_Dofs_3 = Delta1(lim3,:);
Lin_Dofs_1 = Delta(lim1,:);
Lin_Dofs_2 = Delta(lim2,:);
Lin_Dofs_3 = Delta(lim3,:);

R1 = (X1(1) * RWG_Dofs_1) + (X1(2) * Lin_Dofs_1);
R2 = (X2(1) * RWG_Dofs_2) + (X2(2) * Lin_Dofs_2);
R3 = (X3(1) * RWG_Dofs_3) + (X3(2) * Lin_Dofs_3);
% R1(:,1) = 0;
% R2(:,1) = 0;
% R3(:,1) = 0;
% R2(26:50,:) = R2(1:25,:) + R2(26:50,:) + R1(26:end,:);
% R2(51:75,:) = R2(76:end,:) + R2(51:75,:) + R3(26:end,:);
% R2(1:2:end,:) = R2(2:2:end,:) + R2(1:2:end,:);
% R = R1 + R2 + R3;
R1(26:50,:) = R1(26:50,:) + R2(1:25,:);
R3(1:25,:) = R3(1:25,:) + R2(26:50,:);

R = [R1(1:50,:);R3(1:end,:)];
% Ruv(26:2:124,:) = [];
%  [node_coords,quad_Elements] = MeshRectanglularPlate(a(iter),a(iter),8,8) ;
 hold on
%  quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), R(:,1), R(:,2), R(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), X(1)*Delta(:,1), X(1)*Delta(:,2), X(1)*Delta(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), X(2)*Delta1(:,1), X(2)*Delta1(:,2), X(2)*Delta1(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% %  [node_coords,quad_Elements] = MeshRectanglularPlate(a(iter),a(iter),8,8) ;
% %  hold on
% %  quiver3(Ruv(1:26,1),Ruv(1:26,2),Ruv(1:26,3), R(:,1), R(:,2), R(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
%  quiver3(Ruv(lim1,1),Ruv(lim1,2),Ruv(lim1,3), R1(:,1), R1(:,2), R1(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
  quiver3(Ruv(25:2:125,1),Ruv(25:2:125,2),Ruv(25:2:125,3), R2(:,1), R2(:,2), R2(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
%     quiver3(Ruv(lim3,1),Ruv(lim3,2),Ruv(lim3,3), R3(:,1), R3(:,2), R3(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % % quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), Delta2(:,1), Delta2(:,2), Delta2(:,3), 1,'LineWidth', 2,'MaxHeadSize',5);
% % quiver3(Edge3ruv(:,1),Edge3ruv(:,2),Edge3ruv(:,3), Edge1(:,1), Edge1(:,2), Edge1(:,3), 0.5);
% axis equal
% triplot(elem_nodes,node_coords(:,1),node_coords(:,2),'red','LineWidth',linewidth_triplot);
% hold on;
% quiver(coord_xy(:,1),coord_xy(:,2),bf1(:,1),bf1(:,2),'LineWidth',linewidth_quiver);
% axis(axis_limits);
% set(gca,'DataAspectRatio',[1 1 1]); % axes equal;