%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PoWER Large Gas Turbine cycle                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given variables         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Temperatures

T_6 = 1800; % turbine inlet temperature
T_7 = 1200; % recuperator inlet temperature
T_32 = 276; % compressor inlet temperature 
T_10 = 363; % mixing recirculation temperature
T_amb = 288; % ambient temperature

% Pressures

p_amb = 1; % ambient pressure

% Efficiencies/Effectiveness

eta_cp = 0.9; % polytropic compressor efficiency
eta_tp = 0.92; % polytropic turbine efficiency

eff_cup = 0.92; % recuperator effectiveness
eff_cir = 0.9; % recirculator effectiveness

% Pressure drops

pi_hr = 0.97; % pressure drop across recuperator hot side
pi_cr = 0.98; % pressure drop across recuperator cold side

pi_ic = 0.99; % pressure drop across intercooler
pi_ec = 0.99; % pressure drop across evaporator
pi_gc = 0.99; % pressure drop across generator

pi_b = 0.975; % pressure drop across burner

% Gamma values

gam_h = 1.3; % gamma for exhaust gases
gam_m = 1.35; % gamma for mixed gases
gam_c = 1.4; % gamma for cold gases

% Bleed

bl = 0.12; % bleed

% Heat and COP

q_f = 45000000; % heating value of the fuel
cop = 0.7;

% Inputs

pi_21 = 4;
r = 3;

% Fuel-Air ratio

phi = 0.9;

% Cp values

R = 287;
cp_h = (gam_h*R)/(gam_h-1);
cp_m = (gam_m*R)/(gam_m-1);
cp_c = (gam_c*R)/(gam_c-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure track           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_1 = p_amb;

p_2 = pi_21*p_1; % based on input

p_3 = p_2; % no pressure drop across mixer

p_10 = p_2; % no pressure drop across mixer

p_31 = pi_ic*p_3;

p_32 = pi_ec*p_3;

% Now we start calculating all the pressure drops

pi_54 = pi_cr;

pi_65 = pi_b;

pi_76 = (T_7/T_6)^((gam_h)/((gam_h-1)*eta_tp)); % pressure drop across turbine

pi_97 = pi_hr;

pi_1091 = pi_gc;

pi_872 = pi_21; % pressure drop across low pressure turbine

% Going in reverse now

p_91 = p_10/pi_1091; % pressure drop across generator

p_9 = p_91; % no pressure drop across mixer

p_7 = p_9/pi_97; % pressure drop across recuperator hot side

p_6 = p_7/pi_76; % pressure drop across turbine

p_5 = p_6/pi_65; % pressure drop across burner

p_4 = p_5/pi_54; % pressure drop across recuperator cold side

p_72 = p_9; % no pressure drop across mixer

p_8 = p_72/pi_21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature track          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_1 = T_amb;

% polytropic efficiency of low pressure compressor
T_2 = (T_1)*(pi_21^((gam_c-1)/(gam_c*eta_cp))); 

T_3 = (((r*cp_m)/cp_c) - T_10)/(((r*cp_m)/cp_c)-1); % mixing at 10,2 and 3

T_31 = T_3 - (eff_cir*(T_3 - T_amb)); % recirculation effectiveness

% polytropic efficiency of high pressure compressor
T_4 = ((p_4/p_32)^((gam_m-1)/(gam_m*eta_cp)))*T_32; 

% We have a set of three equations to calculate T_9, f_l and T_5
t1 = (gam_h*(gam_m-1))/(gam_m*(gam_h-1));
t3 = (q_f/cp_m)-T_6;
t4 = T_6 - T_4;
s1 = (eff_cup*(T_7 - T_4))/t1;
p1 = T_7*(t3+t4);
p2 = s1*t4 + s1*t1*T_7;
p3 = t3+t4-s1;

T_9 = (p1-p2)/p3;
t2 = T_7 - T_9;

f_l = (t4-t1*t2)/(t3+t1*t2);

T_5 = T_6 - f_l*t3;
% Done with the set of equations

T_91 = T_9; % just a splitting stream
T_72 = T_9;

T_8 = T_72+((1/eta_tp)*(((pi_21)^(gam_m/(gam_m-1)))-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recalculating Low pressure compressor pressure ratio   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi_21_new = ((T_72*eta_cp*eta_tp)/T_1)^(gam_c/(gam_c-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heat and Work          %
%%%%%%%%%%%%%%%%%%%%%%%%%%

q_ge = r*cp_m*(T_10 - T_91); % generator heat
q_ev = (1+r)*cp_m*(T_32 - T_31); % evaporator heat

q_ext = (cop*q_ev)-q_ge; % external load

f_g = (1+r)*f_l; % global fuel-air ratio based on recirculation
q_in = f_g*q_f; % input heat

w = ((0.88*(1+r))+f_l)*cp_h*(T_6-T_7); % work output

eff = (w)/q_in; % total efficiency
