function [] = Format3DGraph(quiverObject,cyl_definition)
% Set the properties of a 3D graph (view angle, quiver arrow colour, etc)
% figure

set(gcf,'Visible','on')
set(gca,'visible','off')

% set(gca,'XTickLabelMode','auto')
% -----------------------%  Thin wire plots % -----------------------------
% semilogx(rho(1:11),imag(G_B_MBF(1:11))*1000,'-s','LineWidth',2,'MarkerSize',8)
% hold on
% semilogx(rho(1:11),imag(G_B(1:11))*1000,':+','LineWidth',2,'MarkerSize',8)
% semilogx(rho(1:11),Tube_Hallen_B(1:11),'-.o','LineWidth',2,'MarkerSize',8)
% semilogx(rho(1:11),Flat_MFIE_B(1:11),'--^','LineWidth',2,'MarkerSize',8)


% -------------------%  Dipole radius label % -----------------------------
% xlabel('Dipole radius, a ($\lambda$)','FontSize',16)
% ylabel('Conductance, G (mS)','FontSize',16)
% legend('MBF, N = 211','RWG without end cap, N = 2020', 'Hallen Equation [REF]', 'MFIE with flat end cap [REF]','Location','northwest','FontSize',14)
% title('Dipole Input Conductance Re$\{1/Z_{in}\}$','FontSize', 20);

% -------------------%  3D View label % -----------------------------------
% xlabel('x ($\lambda$)','FontSize',16)
% ylabel('y ($\lambda$)','FontSize',16)
% zlabel('z ($\lambda$)','FontSize',16)
% set(gca,'FontName','CMU Serif','FontSize',16)






% xlim([0,0.2])
% ylim([0,0.2])
% zlim([0,0.2])
view(-45,30) % Set the view angle
caxis([0 2]) % Set the colourbar limits

% Quiver Properties
quiverObject.Color = 'black';
quiverObject.LineWidth = 1;
quiverObject.AutoScaleFactor = 1;

% Colourbar properties
% a = colorbar;
set(colorbar,'visible','off')
% a.Label.Interpreter = 'latex';
% a.Label.String = '$|\textbf{J}_{MBF}|$';


% ax.Position = [100 100 540 400];

set(gcf,'Position',[100 100 500 500])

hold on 
zAxis = quiver3(0.75,0,0,-1,0,0);
zAxis.Color = 'black';
zAxis.LineWidth = 1;
zAxis.AutoScaleFactor = 1;
textObj = text(-0.35,0,0,"z");

% ax = gca;
textObj.FontSize = 16; 
textObj.FontName = "CMU Serif"; 

if cyl_definition.MBF == "both"
    identifierName = "full_mbf" + "_" + cyl_definition.radius+ cyl_definition.firstNode + "_" + cyl_definition.lastNode + cyl_definition.vertices;
else
   identifierName = cyl_definition.MBF + "_" +  cyl_definition.radius+cyl_definition.firstNode + "_" +cyl_definition.lastNode + cyl_definition.vertices; 
end
path = sprintf("C:\\Users\\jduplessis\\Desktop\\Journal_Article_Figs\\%s.png",identifierName)
print("Exporting figure")
% export_fig path gcf
exportgraphics(gcf,path,'Resolution',500)
end