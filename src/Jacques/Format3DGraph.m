function [] = Format3DGraph()
% Set the properties of a 3D graph (view angle, quiver arrow colour, etc)
% figure
set(gcf,'Visible','on')
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
xlabel('x ($\lambda$)','FontSize',16)
ylabel('y ($\lambda$)','FontSize',16)
zlabel('z ($\lambda$)','FontSize',16)
% set(gca,'FontName','CMU Serif','FontSize',16)

ax = gca;
ax.FontSize = 16; 
ax.FontName = "CMU Serif"; 

% xlim([0,0.2])
% ylim([0,0.2])
% zlim([0,0.2])
view(-45,30) % Set the view angle
caxis([0 2]) % Set the colourbar limits

end