%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wind Tunnel Construction                     								  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Given          %
%%%%%%%%%%%%%%%%%%

T_01 = 300; % Stagnation temperature at inlet
p_amb = 14.7;
t = 180; % time taken for blowdown in seconds

%%%%%%%%%%%%%%%%%
% Given cases   %
%%%%%%%%%%%%%%%%%

h = [4,6,12];
w = [3,6,12];
M = [1.5,2,3,5];
dif_ang = [3,4,5]; % diffuser angles

% Case-wise pressures

p_c = [225,600,1000,2000]; % case pressures

%%%%%%%%%%%%%%%%
% Constants    %
%%%%%%%%%%%%%%%%

gam = 1.4; % gamma
R = 287;
cp = 1005;

%%%%%%%%%%%%%%%
% Conversions %
%%%%%%%%%%%%%%%

in_to_m = 0.0254;
lb_to_kg = 0.453592;
ft_to_m = 0.3058;

%%%%%%%%%%%%%%%%%%%%
% Basic relations  %
%%%%%%%%%%%%%%%%%%%%

A = h.*w;
p_e = p_amb;

%%%%%%%%%%%%%%%%%%
% Part A         %
%%%%%%%%%%%%%%%%%%

% Subset 1

A_t = zeros([3,4]); % throat area

A_by_Astar = (1.2.^(-3)).*(((1+0.2.*(M.^2)).^3)./M); % based on Mach number

for i=1:3
	for j=1:4
		A_t(i,j) = double(A(i))/A_by_Astar(j);
	end
end

% Subset 2

A_s = zeros([3,3]); % shock area
x = zeros([3,3]);
h_s = zeros([3,3]);
w_s = zeros([3,3]);

for i=1:3
	for j=1:3
		x(i,j) = 5*double(h(i))*tand(double(dif_ang(j)));
		h_s(i,j) = double(h(i))+2*x(i,j);
		w_s(i,j) = double(w(i))+2*x(i,j);
		A_s(i,j) = h_s(i,j)*w_s(i,j);
	end
end

As_by_A = zeros([3,3]); % shock area by given area

for i=1:3
	for j=1:3
		As_by_A(i,j) = A_s(i,j)/double(A(i));
	end
end

As_by_Astar = zeros([3,3,4]); % to get Mach number before shock

for i=1:3
	for j=1:3
		for k=1:4
			As_by_Astar(i,j,k) = As_by_A(i,j)*A_by_Astar(k);
		end
	end
end

M_x = zeros(3,3,4); % Mach number upstream of shock
T_x0 = zeros(3,3,4); % static to stagnation temperature ratio
pi_x0 = zeros(3,3,4); % static to stagnation pressure ratio
rho_x0 = zeros(3,3,4); % static to stagnation density ratio
Ax_by_Astar = zeros(3,3,4); % shock area to critical area

for i=1:3
	for j=1:3
		for k=1:4
			[M_x(i,j,k),T_x0(i,j,k),pi_x0(i,j,k),rho_x0(i,j,k),...
				Ax_by_Astar(i,j,k)] = ...
				flowisentropic(1.4,As_by_Astar(i,j,k),'sup');
		end
	end
end

M_up = zeros(3,3,4); % upstream mach number
T_yx = zeros(3,3,4); % static temperature ratio downstream to upstream of shock
pi_yx = zeros(3,3,4); % static pressure ratio downstream to upstream of shock
rho_yx = zeros(3,3,4); % density ratio downstream to upstream of shock
M_y = zeros(3,3,4); % downstream shock mach number
pi_0y0x = zeros (3,3,4); % stagnation pressure ratio downstream to upstream of shock
pi_0yx = zeros(3,3,4); % stagnation pressure downstream to static pressure upstream

for i=1:3
	for j=1:3
		for k=1:4
			[M_up(i,j,k),T_yx(i,j,k),pi_yx(i,j,k),rho_yx(i,j,k),M_y(i,j,k),...
				pi_0y0x(i,j,k),pi_0yx(i,j,k)]...
				= flownormalshock(1.4,M_x(i,j,k),'mach');
		end
	end
end

p_y = zeros(3,3,4); % pressure at the end of shock

for i=1:3
	for j=1:3
		for k=1:4
			p_y(i,j,k) = p_e; % assuming shock is at the exit
		end
	end
end

p_x = p_y./pi_yx;

p_0x = p_x./pi_x0;

p_01 = p_0x; % stagnation pressure in the tank

si_A_t = A_t.*((in_to_m)^2); % in sq.metres
si_p_01 = (p_01./14.7); % in atm
pa_p_01 = si_p_01.*(101324); % in pascal
si_p_c = 101324.*(p_c./14.7); %in pascal

% Subset 3

t1 = zeros(3,3,4);
t2 = zeros(3,3,4);
t3 = zeros(3,3,4);
si_flow_rate = zeros(3,3,4); % mass flow rate of the system

for i=1:3
	for j=1:3
		for k=1:4
			t1(i,j,k) = (pa_p_01(i,j,k))/sqrt(double(R*T_01));
			t2(i,j,k) = double(si_A_t(i,k))*sqrt(gam);
			t3(i,j,k) = (1.2)^(double(gam+1)/double(2-2*gam));
			si_flow_rate(i,j,k) = t1(i,j,k)*t2(i,j,k)*t3(i,j,k);
		end
	end
end

% Subset 4

% Isothermal case

si_V = zeros(3,3,4,4); % volume required
s1 = zeros(3,3,4,4);
s2 = zeros(3,3,4,4);

for i=1:3
	for j=1:3
		for k=1:4
			for l=1:4
				s1 = si_flow_rate(i,j,k)*double(R)*double(t)*double(T_01);
				s2 = si_p_c(l) - pa_p_01(i,j,k);
				si_V(i,j,k,l) = s1/s2;
			end
		end
	end
end

% Isentropic case

T_0c = zeros(3,3,4,4); % stagnation temperature in terms of pressures

for i=1:3
	for j=1:3
		for k=1:4
			for l=1:4
				T_0c(i,j,k,l) = double(T_01)*((pa_p_01(i,j,k)/...
					si_p_c(l))^(double(R)/double(cp)));
			end
		end
	end
end

si_V_s = zeros(3,3,4,4); % isentropic volume required
s1_s = zeros(3,3,4,4);
s2_s = zeros(3,3,4,4);

for i=1:3
	for j=1:3
		for k=1:4
			for l=1:4
				s1_s = si_flow_rate(i,j,k)*double(R)*double(t);
				s2_s = (si_p_c(l)/T_01) - (pa_p_01(i,j,k)/T_0c(i,j,k,l));
				si_V_s(i,j,k,l) = s1_s/s2_s;
			end
		end
	end
end

%%%%%%%%%%%%%%%%%%%%
% Part B           %
%%%%%%%%%%%%%%%%%%%%

% New given conditions

M1 = [1.5,2,3];
a1 = [6,12];
dif_ang1 = 3;

l = 2.*a1;

% Necessary parameters

A1 = a1.^2; % area
P1 = 4.*a1; % perimeter
D_H = (4.*A1)./P1; % hydraulic diameter
epsilon = 0.015e-3; % surface roughness based on the material - stainless steel

% Calculations

pa_p_1 = zeros(2,3); % static pressure at inlet

for i=2:3
	for j=1:3
		p1 = (1 + ((double(gam)-1)/2)*(M1(j)^2))^(double(gam)/(double(gam)-1));
		pa_p_1(i-1,j) = pa_p_01(i,1,j)/p1;
	end
end

T_0t = 1 + ((gam-1)/2).*(double(M1.^2)); % stagnation pressure to static pressure ratio
T_1 = T_01./T_0t; % static pressure in tank
a_sound = sqrt((gam*R).*T_1); % speed of sound
v = a_sound.*M1; % velocity

rho = zeros(2,3); % density based on perfect gas law

for i=1:2
	for j=1:3
		rho(i,j) = pa_p_1(i,j)/(R*T_1(j));
	end
end

mu_0 = 18.27e-6;
T_0_s = 291.15;
C_s = 120;
mu = zeros(1,3); % viscosity based on Sutherland's law

for i=1:3
	mu(i) = mu_0*((double(T_0_s) + double(C_s))/(T_1(i)+...
		double(C_s)))*((T_1(i)/double(T_0_s))^1.5);
end

sr = epsilon./D_H; % surface roughness factor

Re = zeros(2,3); % Reynolds number

for i=1:2
	for j=1:3
		Re(i,j) = (rho(i,j)/mu(j))*v(j)*double(D_H(i));
	end
end

f_D = zeros(2,3); % Darcy friction factor

for i=1:2
	for j=1:3
		f_D(i,j) = 0.0625/...
			(log((sr(i)/3.7)+(5.74/(Re(i,j)^0.9)))^2); % Moody chart equation
	end
end

f = f_D./4; % Fanning friction factor is a quarter of Darcy friction factor

A2 = zeros(2,3); % new area of the duct

for i=1:2
	for j=1:3
		A2(i,j) = A1(i) + A1(i)*((2*double(gam)*f(i,j)*l(i))/double(D_H(i)));
	end
end

a2 = sqrt(A2); % new side length of the square duct

da = zeros(2,3); % change in lengths

for i=1:2
	for j=1:3
		da(i,j) = (a2(i,j)-a1(i));
	end
end

theta = zeros(2,3); % angle to be diverged

for i=1:2
	for j=1:3
		theta(i,j) = atand(da(i,j)/(2*l(i)));
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collecting variables for plots        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stagnation pressure

stag_3_43 = zeros(1,4); % stagnation pressures based on areas and diffuser angles
stag_4_43 = zeros(1,4);
stag_5_43 = zeros(1,4);
stag_3_66 = zeros(1,4);
stag_4_66 = zeros(1,4);
stag_5_66 = zeros(1,4);
stag_3_12 = zeros(1,4);
stag_4_12 = zeros(1,4);
stag_5_12 = zeros(1,4);

for i=1:4
	stag_3_43(i) = pa_p_01(1,1,i);
	stag_4_43(i) = pa_p_01(1,2,i);
	stag_5_43(i) = pa_p_01(1,3,i);
	stag_3_66(i) = pa_p_01(2,1,i);
	stag_4_66(i) = pa_p_01(2,2,i);
	stag_5_66(i) = pa_p_01(2,3,i);
	stag_3_12(i) = pa_p_01(3,1,i);
	stag_4_12(i) = pa_p_01(3,2,i);
	stag_5_12(i) = pa_p_01(3,3,i);
end

% Mass flow rate

flow_3_43 = zeros(1,4); % flow rates collected based on areas and diffuser angles
flow_4_43 = zeros(1,4);
flow_5_43 = zeros(1,4);
flow_3_66 = zeros(1,4);
flow_4_66 = zeros(1,4);
flow_5_66 = zeros(1,4);
flow_3_12 = zeros(1,4);
flow_4_12 = zeros(1,4);
flow_5_12 = zeros(1,4);

for i=1:4
	flow_3_43(i) = si_flow_rate(1,1,i);
	flow_4_43(i) = si_flow_rate(1,2,i);
	flow_5_43(i) = si_flow_rate(1,3,i);
	flow_3_66(i) = si_flow_rate(2,1,i);
	flow_4_66(i) = si_flow_rate(2,2,i);
	flow_5_66(i) = si_flow_rate(2,3,i);
	flow_3_12(i) = si_flow_rate(3,1,i);
	flow_4_12(i) = si_flow_rate(3,2,i);
	flow_5_12(i) = si_flow_rate(3,3,i);
end

% Isothermal blowdown

vol_3_43_15 = zeros(1,4); % case-wise blowdown volumes
vol_4_43_15 = zeros(1,4);
vol_5_43_15 = zeros(1,4);
vol_3_66_15 = zeros(1,4);
vol_4_66_15 = zeros(1,4);
vol_5_66_15 = zeros(1,4);
vol_3_12_15 = zeros(1,4);
vol_4_12_15 = zeros(1,4);
vol_5_12_15 = zeros(1,4);
vol_3_43_2 = zeros(1,4);
vol_4_43_2 = zeros(1,4);
vol_5_43_2 = zeros(1,4);
vol_3_66_2 = zeros(1,4);
vol_4_66_2 = zeros(1,4);
vol_5_66_2 = zeros(1,4);
vol_3_12_2 = zeros(1,4);
vol_4_12_2 = zeros(1,4);
vol_5_12_2 = zeros(1,4);
vol_3_43_3 = zeros(1,4);
vol_4_43_3 = zeros(1,4);
vol_5_43_3 = zeros(1,4);
vol_3_66_3 = zeros(1,4);
vol_4_66_3 = zeros(1,4);
vol_5_66_3 = zeros(1,4);
vol_3_12_3 = zeros(1,4);
vol_4_12_3 = zeros(1,4);
vol_5_12_3 = zeros(1,4);

for i=1:4
	vol_3_43_15(i) = si_V(1,1,1,i);
	vol_4_43_15(i) = si_V(1,2,1,i);
	vol_5_43_15(i) = si_V(1,3,1,i);
	vol_3_66_15(i) = si_V(2,1,1,i);
	vol_4_66_15(i) = si_V(2,2,1,i);
	vol_5_66_15(i) = si_V(2,3,1,i);
	vol_3_12_15(i) = si_V(3,1,1,i);
	vol_4_12_15(i) = si_V(3,2,1,i);
	vol_5_12_15(i) = si_V(3,3,1,i);
	vol_3_43_2(i) = si_V(1,1,2,i);
	vol_4_43_2(i) = si_V(1,2,2,i);
	vol_5_43_2(i) = si_V(1,3,2,i);
	vol_3_66_2(i) = si_V(2,1,2,i);
	vol_4_66_2(i) = si_V(2,2,2,i);
	vol_5_66_2(i) = si_V(2,3,2,i);
	vol_3_12_2(i) = si_V(3,1,2,i);
	vol_4_12_2(i) = si_V(3,2,2,i);
	vol_5_12_2(i) = si_V(3,3,2,i);
	vol_3_43_3(i) = si_V(1,1,3,i);
	vol_4_43_3(i) = si_V(1,2,3,i);
	vol_5_43_3(i) = si_V(1,3,3,i);
	vol_3_66_3(i) = si_V(2,1,3,i);
	vol_4_66_3(i) = si_V(2,2,3,i);
	vol_5_66_3(i) = si_V(2,3,3,i);
	vol_3_12_3(i) = si_V(3,1,3,i);
	vol_4_12_3(i) = si_V(3,2,3,i);
	vol_5_12_3(i) = si_V(3,3,3,i);
end

% Isothermal blowdown for Mach 5

vol_3_43_5 = zeros(1,4); % case-wise blowdown volumes
vol_4_43_5 = zeros(1,4);
vol_5_43_5 = zeros(1,4);
vol_3_66_5 = zeros(1,4);
vol_4_66_5 = zeros(1,4);
vol_5_66_5 = zeros(1,4);
vol_3_12_5 = zeros(1,4);
vol_4_12_5 = zeros(1,4);
vol_5_12_5 = zeros(1,4);

for i=3:4
	vol_3_43_5(i) = si_V(1,1,4,i);
	vol_4_43_5(i) = si_V(1,2,4,i);
	vol_5_43_5(i) = si_V(1,3,4,i);
	vol_3_66_5(i) = si_V(2,1,4,i);
	vol_4_66_5(i) = si_V(2,2,4,i);
	vol_5_66_5(i) = si_V(2,3,4,i);
	vol_3_12_5(i) = si_V(3,1,4,i);
	vol_4_12_5(i) = si_V(3,2,4,i);
	vol_5_12_5(i) = si_V(3,3,4,i);
end

% Isentropic blowdown

svol_3_43_15 = zeros(1,4); % case-wise blowdown volumes
svol_4_43_15 = zeros(1,4);
svol_5_43_15 = zeros(1,4);
svol_3_66_15 = zeros(1,4);
svol_4_66_15 = zeros(1,4);
svol_5_66_15 = zeros(1,4);
svol_3_12_15 = zeros(1,4);
svol_4_12_15 = zeros(1,4);
svol_5_12_15 = zeros(1,4);
svol_3_43_2 = zeros(1,4);
svol_4_43_2 = zeros(1,4);
svol_5_43_2 = zeros(1,4);
svol_3_66_2 = zeros(1,4);
svol_4_66_2 = zeros(1,4);
svol_5_66_2 = zeros(1,4);
svol_3_12_2 = zeros(1,4);
svol_4_12_2 = zeros(1,4);
svol_5_12_2 = zeros(1,4);
svol_3_43_3 = zeros(1,4);
svol_4_43_3 = zeros(1,4);
svol_5_43_3 = zeros(1,4);
svol_3_66_3 = zeros(1,4);
svol_4_66_3 = zeros(1,4);
svol_5_66_3 = zeros(1,4);
svol_3_12_3 = zeros(1,4);
svol_4_12_3 = zeros(1,4);
svol_5_12_3 = zeros(1,4);

for i=1:4
	svol_3_43_15(i) = si_V_s(1,1,1,i);
	svol_4_43_15(i) = si_V_s(1,2,1,i);
	svol_5_43_15(i) = si_V_s(1,3,1,i);
	svol_3_66_15(i) = si_V_s(2,1,1,i);
	svol_4_66_15(i) = si_V_s(2,2,1,i);
	svol_5_66_15(i) = si_V_s(2,3,1,i);
	svol_3_12_15(i) = si_V_s(3,1,1,i);
	svol_4_12_15(i) = si_V_s(3,2,1,i);
	svol_5_12_15(i) = si_V_s(3,3,1,i);
	svol_3_43_2(i) = si_V_s(1,1,2,i);
	svol_4_43_2(i) = si_V_s(1,2,2,i);
	svol_5_43_2(i) = si_V_s(1,3,2,i);
	svol_3_66_2(i) = si_V_s(2,1,2,i);
	svol_4_66_2(i) = si_V_s(2,2,2,i);
	svol_5_66_2(i) = si_V_s(2,3,2,i);
	svol_3_12_2(i) = si_V_s(3,1,2,i);
	svol_4_12_2(i) = si_V_s(3,2,2,i);
	svol_5_12_2(i) = si_V_s(3,3,2,i);
	svol_3_43_3(i) = si_V_s(1,1,3,i);
	svol_4_43_3(i) = si_V_s(1,2,3,i);
	svol_5_43_3(i) = si_V_s(1,3,3,i);
	svol_3_66_3(i) = si_V_s(2,1,3,i);
	svol_4_66_3(i) = si_V_s(2,2,3,i);
	svol_5_66_3(i) = si_V_s(2,3,3,i);
	svol_3_12_3(i) = si_V_s(3,1,3,i);
	svol_4_12_3(i) = si_V_s(3,2,3,i);
	svol_5_12_3(i) = si_V_s(3,3,3,i);
end

% Isentropic blowdown for Mach 5

svol_3_43_5 = zeros(1,4); % case-wise isentropic blowdown volumes
svol_4_43_5 = zeros(1,4);
svol_5_43_5 = zeros(1,4);
svol_3_66_5 = zeros(1,4);
svol_4_66_5 = zeros(1,4);
svol_5_66_5 = zeros(1,4);
svol_3_12_5 = zeros(1,4);
svol_4_12_5 = zeros(1,4);
svol_5_12_5 = zeros(1,4);

for i=3:4
	svol_3_43_5(i) = si_V_s(1,1,4,i);
	svol_4_43_5(i) = si_V_s(1,2,4,i);
	svol_5_43_5(i) = si_V_s(1,3,4,i);
	svol_3_66_5(i) = si_V_s(2,1,4,i);
	svol_4_66_5(i) = si_V_s(2,2,4,i);
	svol_5_66_5(i) = si_V_s(2,3,4,i);
	svol_3_12_5(i) = si_V_s(3,1,4,i);
	svol_4_12_5(i) = si_V_s(3,2,4,i);
	svol_5_12_5(i) = si_V_s(3,3,4,i);
end

%%%%%%%%%%%%%%%%%%;
% Plots           %
%%%%%%%%%%%%%%%%%%%

% Throat area

figure('visible','off')
plot1 = plot(M,A_t(1,:),'g-*',M,A_t(2,:),'b-*',M,A_t(3,:),'r-*');
legend('A = 4x3','A = 6x6','A = 12x12');
xlabel('Mach number');
ylabel('Throat area (in sq.inches)');
title('Throat area vs Mach number');
print('-depsc','./images/throat_area_vs_mach');

% Stagnation pressure

figure('visible','off')
plot2 = plot(M,stag_3_43,'-*',M,stag_4_43,'-*',M,stag_5_43,'-*',...
	M,stag_3_66,'--*',M,stag_4_66,'--*',M,stag_5_66,'--*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'Location','northwest','boxoff');
xlabel('Mach number');
ylabel('Stagnation pressure (in Pa)');
title('Stagnation pressure vs Mach number');
print('-depsc','./images/stagnation_pressure_vs_mach');

figure('visible','off')
plot2_1 = plot(M,stag_3_12,'-*',M,stag_4_12,'-*',M,stag_5_12,'-*');
legend('12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northwest','boxoff');
xlabel('Mach number');
ylabel('Stagnation pressure for 12x12 cross section (in Pa)');
title('Stagnation pressure vs Mach number');
print('-depsc','./images/stag_12x12_vs_mach');

% Mass flow rate

figure('visible','off')
plot_3 = plot(M,flow_3_43,'-*',M,flow_4_43,'-*',M,flow_5_43,'-*',...
	M,flow_3_66,'--*',M,flow_4_66,'--*',M,flow_5_66,'--*',...
	M,flow_3_12,'-.*',M,flow_4_12,'-.*',M,flow_5_12,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','eastoutside','boxoff');
xlabel('Mach number');
ylabel('Mass flow rate (in kg/s)');
title('Mass flow rate vs Mach number');
print('-depsc','./images/flow_vs_mach');

% Isothermal blowdown

figure('visible','off')
plot4_1 = plot(si_p_c,vol_3_43_15,'-*',si_p_c,vol_4_43_15,'-*',...
	si_p_c,vol_5_43_15,'-*',...
	si_p_c,vol_3_66_15,'--*',si_p_c,vol_4_66_15,'--*',si_p_c,vol_5_66_15,'--*',...
	si_p_c,vol_3_12_15,'-.*',si_p_c,vol_4_12_15,'-.*',si_p_c,vol_5_12_15,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 1.5');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/vol_vs_case_mach_15');

figure('visible','off')
plot4_2 = plot(si_p_c,vol_3_43_2,'-*',si_p_c,vol_4_43_2,'-*',...
	si_p_c,vol_5_43_2,'-*',...
	si_p_c,vol_3_66_2,'--*',si_p_c,vol_4_66_2,'--*',si_p_c,vol_5_66_2,'--*',...
	si_p_c,vol_3_12_2,'-.*',si_p_c,vol_4_12_2,'-.*',si_p_c,vol_5_12_2,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 2');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/vol_vs_case_mach_2');

figure('visible','off')
plot4_3 = plot(si_p_c,vol_3_43_3,'-*',si_p_c,vol_4_43_3,'-*',...
	si_p_c,vol_5_43_3,'-*',...
	si_p_c,vol_3_66_3,'--*',si_p_c,vol_4_66_3,'--*',si_p_c,vol_5_66_3,'--*',...
	si_p_c,vol_3_12_3,'-.*',si_p_c,vol_4_12_3,'-.*',si_p_c,vol_5_12_3,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 3');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/vol_vs_case_mach_3');

figure('visible','off')
plot6_2 = plot(si_p_c(3:4),vol_3_43_5(3:4),'-*',si_p_c(3:4),vol_4_43_3(3:4),'-*',...
	si_p_c(3:4),vol_5_43_3(3:4),'-*',...
	si_p_c(3:4),vol_3_66_3(3:4),'--*',si_p_c(3:4),vol_4_66_3(3:4),'--*',...
	si_p_c(3:4),vol_5_66_3(3:4),'--*',...
	si_p_c(3:4),vol_3_12_3(3:4),'-.*',si_p_c(3:4),vol_4_12_3(3:4),'-.*',...
	si_p_c(3:4),vol_5_12_3(3:4),'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 5');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/vol_vs_case_mach_5');

% Isentropic blowdown

figure('visible','off')
plot5_1 = plot(si_p_c,svol_3_43_15,'-*',si_p_c,svol_4_43_15,'-*',...
	si_p_c,svol_5_43_15,'-*',...
	si_p_c,svol_3_66_15,'--*',si_p_c,svol_4_66_15,'--*',si_p_c,svol_5_66_15,'--*',...
	si_p_c,svol_3_12_15,'-.*',si_p_c,svol_4_12_15,'-.*',si_p_c,svol_5_12_15,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 1.5');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/svol_vs_case_mach_15');

figure('visible','off')
plot5_2 = plot(si_p_c,svol_3_43_2,'-*',si_p_c,svol_4_43_2,'-*',...
	si_p_c,svol_5_43_2,'-*',...
	si_p_c,svol_3_66_2,'--*',si_p_c,svol_4_66_2,'--*',si_p_c,svol_5_66_2,'--*',...
	si_p_c,svol_3_12_2,'-.*',si_p_c,svol_4_12_2,'-.*',si_p_c,svol_5_12_2,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 2');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/svol_vs_case_mach_2');

figure('visible','off')
plot5_3 = plot(si_p_c,svol_3_43_3,'-*',si_p_c,svol_4_43_3,'-*',...
	si_p_c,svol_5_43_3,'-*',...
	si_p_c,svol_3_66_3,'--*',si_p_c,svol_4_66_3,'--*',si_p_c,svol_5_66_3,'--*',...
	si_p_c,svol_3_12_3,'-.*',si_p_c,svol_4_12_3,'-.*',si_p_c,svol_5_12_3,'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 3');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/svol_vs_case_mach_3');

figure('visible','off')
plot6_1 = plot(si_p_c(3:4),svol_3_43_5(3:4),'-*',...
	si_p_c(3:4),svol_4_43_3(3:4),'-*',...
	si_p_c(3:4),svol_5_43_3(3:4),'-*',...
	si_p_c(3:4),svol_3_66_3(3:4),'--*',si_p_c(3:4),svol_4_66_3(3:4),'--*',...
	si_p_c(3:4),svol_5_66_3(3:4),'--*',...
	si_p_c(3:4),svol_3_12_3(3:4),'-.*',si_p_c(3:4),svol_4_12_3(3:4),'-.*',...
	si_p_c(3:4),svol_5_12_3(3:4),'-.*');
legend('4x3 and 3 degrees','4x3 and 4 degrees','4x3 and 5 degrees',...
	'6x6 and 3 degrees','6x6 and 4 degrees','6x6 and 5 degrees',...
	'12x12 and 3 degrees','12x12 and 4 degrees','12x12 and 5 degrees',...
	'Location','northeast','boxoff');
xlabel('Case Pressures (in Pa), M = 5');
ylabel('Volume required (in cubic metres)');
title('Case pressures vs Volume required for a given Mach number');
print('-depsc','./images/svol_vs_case_mach_5');

% Part B

figure('visible','off')
plot7 = plot(M1,theta(1,:),'-*');
legend('A = 6x6');
xlabel('Mach number, Area = 6x6');
ylabel('Angle to be deviated by');
title('Angle of test section vs Mach number');
print('-depsc','./images/angle_66_vs_mach');

figure('visible','off')
plot8 = plot(M1,theta(2,:),'-*');
legend('A = 12x12');
xlabel('Mach number, Area = 12x12');
ylabel('Angle to be deviated by');
title('Angle of test section vs Mach number');
print('-depsc','./images/angle_12_vs_mach');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The End                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
