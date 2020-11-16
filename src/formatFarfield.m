function formatFarfield(ff_mbf, ff_ref, n_mbf, n_ref)


% Create figure
figure('OuterPosition',[105 233 1630 769]);

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create multiple lines using matrix input to semilogy
% semilogy1 = semilogy(X1,YMatrix1,'LineWidth',1.5);
semilogy_ref = semilogy(ff_ref(:,1),sqrt(ff_ref(:,3).^2 + ff_ref(:,5).^2),'LineWidth',1.5);
semilogy_mbf = semilogy(ff_mbf(:,1),sqrt(ff_mbf(:,3).^2 + ff_mbf(:,5).^2),'LineWidth',1.5);

str_ff1 = sprintf("MBF, N = %i",n_mbf);
str_ff2 = sprintf("RWG, N = %i",n_ref);
set(semilogy_mbf,'DisplayName',str_ff1,'MarkerSize',4,'Marker','o',...
    'Color',[0 0.447058823529412 0.741176470588235]);
set(semilogy_ref,'DisplayName',str_ff2,'Color',[0 0 0]);


% Create ylabel
ylabel({'$|\mathbf{E}_\theta|$ (V)'});

% Create xlabel
xlabel({'$\theta$ (rad)'});

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 6.28318530717959]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','CMU Serif','FontSize',20,'YMinorTick','on','YScale',...
    'log');
% Create legend
legend(axes1,'show');
grid on;
xlim([0 2*pi]);

