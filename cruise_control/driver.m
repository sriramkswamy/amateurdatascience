clear all;
clc;

v_star = 30; % equilibrium value

% Given parameters
m = 1500;
beta_drag = 0.324;
tau_1 = 50;
tau_2 = 5;
k = 4774.65;

% From the derived state space model
A_plant = [-(2*beta_drag*v_star)/m 1/m 0;
			0 0 1;
			0 -1/(tau_1*tau_2) -((1/tau_1)+(1/tau_2));];

B_plant = [0;0;k/(tau_1*tau_2)];

C_plant = [1 0 0];

D_plant = [0];

[num_plant, den_plant] = ss2tf(A_plant,B_plant,C_plant,D_plant);

P = tf(num_plant,den_plant); % The resulting linearized plant model

% Performance specifications given
w_r = 0.01;
eps_r = 0.01;
w_n = 20;
eps_n = 0.01;
w_gc_des = 0.2;
w_pm_des = 30;

% The following code upto the point specified was based on the original loop shaping example discussed in class

% Define P_lin(jw)
Pjw = @(w)(0.1273/((5.184e-5 - 0.233*(w^2))+(sqrt(-1)*(0.006851*w - (w^3)))));

% Make an initial poor P controller based on the Performance parameters
k_p1 = (1/eps_r+1)/abs(Pjw(w_r));
k_p2 = 1/abs(Pjw(w_gc_des));
k_p = min(k_p2,k_p1);
L1 = k_p*P;
C1 = k_p;

% Define the logspace we intend to operate
omegas = logspace(-3,3,100);

% Plot the resulting P controller - useful while designing
% loopfigure = figure;
% bodemag(L1,omegas)
% grid on
% title('P-controller')

% We notice that P controller doesn't really satisfy any condition except noise rejection. So, we try to add a lag controller
% Parameters based on trial and error
p_lag = 10;
z_lag = 1000;
klag = 101/abs(Pjw(0)*k_p); % Make an initial guess for lag gain based on the fact that this will satisfy the given criterion if it satisfies (given criterion + 1)
C_lag = zpk(-z_lag,-p_lag,15*klag*p_lag/z_lag);

% Find the minimal realization of the resulting transfer function
L2 = minreal(series(L1,C_lag));
C2 = series(C1,C_lag);

% Plot the P controller and Lag controller side-by-side - useful while designing
% figure(loopfigure); hold on
% bodemag(L2,omegas,'r--')
% grid on
% set(gca,'ytick',[-100:20:40])
% legend('P-controller','P + Lag controller')

% We see that the controller still doesn't satisfy our requirements. The gain crossover frequency and phase margins are unsatisfied.
% We now add a lead controller to the previous controller to achieve this.
% Parameters based on trial and error
p_lead = 5;
z_lead = 0.1;
C_lead = zpk(-z_lead,-p_lead,1*p_lead/z_lead);

% L3 = minreal(series(series(L2,C_lead),1)); % A single controller doesn't satisfy the phase margin performance
L3 = minreal(series(series(L2,C_lead),C_lead)); % Using of two lead controllers does satisfy
C3 = series(C2,C_lead);

% Plot the resulting controller - useful while designing
% figure
% bodemag(L3,omegas)
% grid on
% set(gca,'ytick',[-100:20:40])
% title('P + lag + 2*lead controller')

% Verify the phase and gain margins - useful while designing
% figure
% margin(L3)

C = minreal(C3);
% Bode plot of the resulting controller - useful while designing
% figure
% bode(C)
% grid on;

% The code upto this part was more or less based on the original code by Dr. Barooah which was used to discuss a loop shaping example in class

% Define the parameters necessary for the plots
t = [0:1:600];

% Define reference based on the control loop (r - y_eq)
r=zeros(length(t),1);
for i=1:1:length(t)
	if t(i)<300
		r(i) = 1.111;
	else
		r(i) = 3.333;
	end
end

% The reference that is given to us
v_des=zeros(length(t),1);
for i=1:1:length(t)
	if t(i)<300
		v_des(i) = 31.111;
	else
		v_des(i) = 33.333;
	end
end

% Find the transfer function for output y and the output
Hyr = minreal((C*P)/(1 + C*P));
yr = lsim(Hyr,r,t);

% Actual output = output + y_eq
v = 30 + yr;

% Find the transfer function for input u and the input
Hur = minreal(C/(1 + C*P));
ur = lsim(Hur,r,t);

% Plot the output y and the reference r
figure
plot(t,yr,'r--',t,r,'b-');
xlabel('Time(s)');
ylabel('Magnitude');
legend('Output - y','Reference - r');
title('Output and reference vs time');

% Plot the input u based on the reference r
figure
plot(t,ur);
xlabel('Time(s)');
ylabel('Magnitude');
legend('Input - u');
title('Input vs time');

% Plot the actual velocity and the desired velocity
figure
plot(t,v,'r--',t,v_des,'b-');
xlabel('Time(s)');
ylabel('Velocity Magnitude');
legend('Actual velocity','Desired velocity');
title('Actual velocity and Desired velocity plot vs time');
