rho = 0.5;
vertices = 12; % Number of vertices
%  Contour = [0 0 0; 1 0 0; 2 0 0];
%    Contour = [0 0 0; 0.5 0 0; 1 0 0; 1.5 0 0; 2 0 0];
%    Contour = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 0.1; 5 0 0.2; 6 0.1 0.3; 7 0.2 0.4; 8 0.25 0.5; 9 0.2 0.6; 10 0.2 0.5; 11 0.2 0.4 ]; % Hardcode some contour
 Contour = [0 0 0; 1 0 0; 2 0 0];

 Contour = RefineMesh(Contour,1);
[points1,quadElements1, elements1] = QuadMesh_v5(Contour,vertices,rho, 0);
[revised_triangles1,triangles1] = QuadtoTri(elements1,vertices, 0);
  PlotMesh(points1,elements1,triangles1, 1);

%   Delta = importDelta('C:\Users\19083688\Desktop\Delta.txt', 1, 1500);
%Delta = importDelta('C:\Users\19083688\Desktop\MoM Codes\V_3.6_feat_plane_wave_MBF_qmom\Delta.txt', 1, 1500);
% ruv = importruv('C:\Users\19083688\Desktop\ruv.txt', 1, 1500);
 Delta = importDeltaNew('C:\Users\19083688\Desktop\Masters\MoMLib\Delta.txt',[1, Inf]);
 Ruv = importruv("C:\Users\19083688\Desktop\Masters\MoMLib\Ruv.txt", [1, Inf]);
 Delta1 = importDeltaNew('C:\Users\19083688\Desktop\Masters\MoMLib\Delta1.txt',[1, Inf]);
 Ruv1 = importruv("C:\Users\19083688\Desktop\Masters\MoMLib\Ruv1.txt", [1, Inf]);
 Delta2 = Delta1 + Delta;
%figure;
%h = subplot(1,3,1);
%   hold on
% quiver3(ruv_large(:,1),ruv_large(:,2),ruv_large(:,3), DeltaLarge(:,1), DeltaLarge(:,2), DeltaLarge(:,3), 5);
% quiver3(Ruv1(:,1),Ruv1(:,2),Ruv1(:,3), Delta1(:,1), Delta1(:,2), Delta1(:,3), 1.5,'LineWidth', 1.5,'MaxHeadSize',3);
% figure
quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), Delta(:,1), Delta(:,2), Delta(:,3), 1,'LineWidth', 1.5,'MaxHeadSize',3);
% quiver3(Ruv(:,1),Ruv(:,2),Ruv(:,3), Delta2(:,1), Delta2(:,2), Delta2(:,3), 1,'LineWidth', 1.5,'MaxHeadSize',3);
% quiver3(Edge3ruv(:,1),Edge3ruv(:,2),Edge3ruv(:,3), Edge1(:,1), Edge1(:,2), Edge1(:,3), 0.5);
axis equal
% triplot(elem_nodes,node_coords(:,1),node_coords(:,2),'red','LineWidth',linewidth_triplot);
% hold on;
% quiver(coord_xy(:,1),coord_xy(:,2),bf1(:,1),bf1(:,2),'LineWidth',linewidth_quiver);
% axis(axis_limits);
% set(gca,'DataAspectRatio',[1 1 1]); % axes equal;