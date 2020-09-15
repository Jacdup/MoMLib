
z = linspace(0,1,20);
deg = linspace(0,2*pi,20);
MBF = sin(deg);
zn = 0.5;
del = z(2) - z(1);
rooftop = ((del-abs(z-zn))/del);
rooftop = rooftop + max(abs(rooftop));
rooftop = rooftop/max(abs(rooftop));
In_axial = kron(MBF,rooftop');
In_axial = reshape(In_axial,[20,20]);
% In_circ  = kron(MBF,MBF);
% In_circ = reshape(In_circ,[20,20]);
In_circ = repmat(MBF,20);
In_circ = In_circ(1,:);
In_circ = reshape(In_circ,[20,20]);
surf(deg, z, In_axial)
figure
surf(deg, z, In_circ');

% Both
figure
% surf(deg,z,In_axial+In_circ')
pcolor(deg,z,In_axial+In_circ')
hold on
h = quiver(x,y,In_circ',In_axial,0.2);



% Azimuthal
figure
ah = axes;
pcolor(deg,z*scale_factor,In_circ')
hold on
h = quiver(x,y*scale_factor,In_circ'*scale_factor,zeros(20),0.2);
xlim([min(x) max(x)])
ylim([min(y) max(y)]*scale_factor)
set(ah,'YTick',y*scale_factor)
set(ah,'YTickLabel',y)

% Axial
u = repmat(rooftop,20);
u =  u(:,1:20);
u = zeros(20);

% v = repmat(v,20);
% v = v(:,1:20);
w = zeros(20);
%  quiver3(x,y,In,u,v,w)
scale_factor = range(x)/range(y);

figure;
ah = axes;
pcolor(deg,z*scale_factor,In_axial)
hold on
h = quiver(x,y*scale_factor,u,v'*scale_factor,0.2);
xlim([min(x) max(x)])
ylim([min(y) max(y)]*scale_factor)
set(ah,'YTick',y*scale_factor)
set(ah,'YTickLabel',y)
% hold on


h.Head.LineStyle = 'solid';