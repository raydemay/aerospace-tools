%% J2 Perturbation Test
clear
clc

%Initialize
mu = 3.986004e14;
R_Earth = 6378.137e3;
%orbital params in inertial frame
r0 = [-2.491984247928895e+06;4.793000892455519e+05;-6.824701828788767e+06];
v0 = [6.708662359765611e+03;2.089444378549832e+03;-2.302871311065931e+03];
f_inc = linspace(0,2*pi);
J2 = 1.08263e-3; %2nd zonal harmonic
%times = [0, 1, 2, 4, 6, 10, 15, 25, 50, 75, 100, 150];
times = linspace(0,365.25*24);
times_sec = times*3600;
rI = zeros(3,length(f_inc));
r_OF = [];

%Get initial orbit conditions
orbit_params = IFtoOF(mu,r0,v0);
a_0 = orbit_params(1);
e_0 = orbit_params(2);
i_0 = orbit_params(3);
LAN_0 = orbit_params(4);
omega_0 = orbit_params(5);
f_0 = orbit_params(6);
p_0 = a_0*(1-e_0^2);
r1 = p_0./(1+e_0*cos(f_inc));
ro = [r1.*cos(f_inc); r1.*sin(f_inc); zeros(1,length(f_inc))];
ROI1 = ROI(omega_0,i_0,LAN_0);
RIO = ROI1';
for j = 1:length(f_inc)
    rI(:,j) = (RIO*ro(:,j))';
end

%Figuring out avg orbit propagation from initial values
T_0 = 2*pi*sqrt(a_0^3/mu);
T_hr = T_0/3600;
K2_0 = (-1.5)*(sqrt(mu)*J2*R_Earth^2)/((1-e_0^2)^2*a_0^(7/2));
omega_dot_avg = K2_0*(2.5*sin(i_0)^2-2);
LAN_dot_avg = K2_0*cos(i_0);
LAN_dot = LAN_dot_avg*T_0;

%Numerically integrate
init_vals = [r0 v0];
tols = odeset('RelTol',1e-12,'AbsTol',1e-15);
[t,x] = ode89(@orbiteqJ2_Earth,times_sec,init_vals,tols);

for k=1:length(times_sec)
    orbit(k,:) = IFtoOF(mu,x(k,1:3),x(k,4:6));
end
p = orbit(:,1).*(1-orbit(:,2).^2);
r = p./(1+orbit(:,2)*cos(f_inc));
r_cell = mat2cell(r,ones(1,length(times_sec)));

for l=1:length(r_cell)
    r_OF{l} = [r_cell{l}.*cos(f_inc); r_cell{l}.*sin(f_inc);zeros(1,length(f_inc))];
end

for m =1:length(orbit(:,1))
    B{m} = ROI(orbit(m,5),orbit(m,3),orbit(m,4));
    rots = cellfun(@transpose,B,'UniformOutput',false);
end

for n=1:length(r_OF)
    for o=1:length(f_inc)
        r_IF{n}(o,:) = rots{n}*r_OF{n}(:,o);
    end
end

%Make an Earth to plot
[c,y,z] = sphere;
c = c*R_Earth;
y = y*R_Earth;
z = z*R_Earth;

figure(1)
plot(times,orbit(:,4))
hold on
title('LAN Variance')
xlabel('Time (Hours)')
ylabel('Longitude of ascending node (radians)')

fig = figure(2);
%fig.WindowState = 'maximized';
plot3(0, 0, 0 , 'bo', 'Markersize', 30, 'Markerfacecolor', 'b');
hold on
surf(c,y,z,'FaceColor','b')
xlim([-8.5e6 8.5e6])
ylim([-8.5e6 8.5e6])
zlim([-8.5e6 8.5e6])
hold on
for q=1:length(r_IF)
    plot3(r_IF{q}(:,1),r_IF{q}(:,2),r_IF{q}(:,3),'LineWidth',2,'Color','r')
    pause(.2)
    drawnow
end

%%  Functions
function derivs = orbiteqJ2_Earth(~,input)
    mu = 3.986004e14;
    R_Earth = 6378.137e3;
    J2 = 1.08263e-3; %2nd zonal harmonic
    derivs = zeros(6,1);
    r = norm(input(1:3));
    derivs(1:3) = input(4:6);
    derivs(4) = -mu*input(1)/r^3+((1.5)*(J2*mu*R_Earth^2/r^4)*(input(1)/r)*(5*(input(3)^2/r^2)-1));
    derivs(5) = -mu*input(2)/r^3+((1.5)*(J2*mu*R_Earth^2/r^4)*(input(2)/r)*(5*(input(3)^2/r^2)-1));
    derivs(6) = -mu*input(3)/r^3+((1.5)*(J2*mu*R_Earth^2/r^4)*(input(3)/r)*(5*(input(3)^2/r^2)-3));
end

function OrbitElems = IFtoOF(mu, r, v)
%Calculate orbital parameters from given inputs
%Get initial unit vectors and "easy" parameters
h = cross (r,v); %angular momentum vector
i_h = h/norm(h); %unit vector for h
c = cross(v,h) - (mu/norm(r))*r; %eccentricity vector
e = norm(c)/mu; %eccentricity
i_e = c/norm(c);
i_y = cross(i_h,i_e);
energy = (norm(v)^2)/2 - mu/norm(r);
a = -mu/(2*energy); %semimajor axis
i = acos(i_h(3)); %inclination in degrees
%Longitude of Ascending node
cos_LAN = -i_h(2)/sin(i);
sin_LAN = i_h(1)/sin(i);
LAN = atan2(sin_LAN, cos_LAN);
%Argument of periapsis
cos_omega = i_y(3)/sin(i);
sin_omega = i_e(3)/sin(i);
omega = atan2(sin_omega, cos_omega);
%DCM from I to O
Roi = [i_e';i_y';i_h']; %single quote is MATLAB transpose
%True anomaly
r_o = Roi*r;
cos_f = r_o(1)/norm(r_o);
sin_f = r_o(2)/norm(r_o);
f = atan2(sin_f, cos_f);
OrbitElems = [a;e;i;LAN;omega;f];
end

function rotations = ROI(om, i, OM)
ie1 = cos(om)*cos(OM)-sin(om)*cos(i)*sin(OM);
ie2 = cos(om)*sin(OM)+sin(om)*cos(i)*cos(OM);
ie3 = sin(om)*sin(i);
iy1 = -(sin(om)*cos(OM)+cos(om)*cos(i)*sin(OM));
iy2 = -sin(om)*sin(OM)+cos(om)*cos(i)*cos(OM);
iy3 = cos(om)*sin(i);
ih1 = sin(i)*sin(OM);
ih2 = -sin(i)*cos(OM);
ih3 = cos(i);
rotations = [ie1 ie2 ie3; iy1 iy2 iy3; ih1 ih2 ih3];
end
