#############################################################################
# PoWER Gas Turbines cycle                                                  #
#############################################################################

#################
# Modules       #
#################

import numpy as np
import matplotlib.pyplot as plt

###################
# Given           #
###################

########## Large Engine ###########

# Temperatures

T6_ld = 1800.0 # HPT inlet

# Polytropic Efficiencies

eta_cp_ld = 0.9 # HPC/LPC 
eta_tp_ld = 0.92 # HPT/LPT 

# HX Effectiveness

eff_cup_ld = 0.92 # RC 
eff_cir_ld = 0.9 # IC/EVA/GEN 

# Pressure drops

pi_cir_ld = 0.99 # IC/EVA/GEN

pi_b_ld = 0.975 # Burner

# Bleed

bl = 0.12

########## Common parameters #########

# Ambient conditions

T_amb = 288.0
p_amb = 1.0

# Fuel-air ratio

phi = 0.9
fg = 0.0611

# Constants

gam_h = 1.3
gam_m = 1.35
gam_c = 1.4
R = 287.0

# Cp values

cp_h = (gam_h*R)/(gam_h-1)
cp_m = (gam_m*R)/(gam_m-1)
cp_c = (gam_c*R)/(gam_c-1)

# Pressure drops

pi_hrd = 0.97
pi_crd = 0.98

# Temperatures

T7d = 1200.0 # RC 
T32d = 276.0 # HPC inlet
T10d = 363.0 # GEN outlet

# Heating value of the fuel

q_f = 45000000.0

# VARS COP

cop = 0.7

########### Small Engine ###########

# Temperatures

T6_sd = 1400.0 # HPT inlet

# Polytropic Efficiencies

eta_cp_sd = 0.82 # HPC/LPC
eta_tp_sd = 0.85 # HPT/LPT

# HX Effectiveness

eff_cup_sd = 0.9 # RC
eff_cir_sd = 0.88 # IC/EVA/GEN

# Pressure drops

pi_b_sd = 0.965 # Burner
pi_cir_sd = 0.985 # IC/EVA/GEN

#######################
# Functions           #
#######################

# Pressure ratio based on polytropic efficiency and Temperatures

def polypic(T1,T2,eta,gam): 
    """
    T1: inlet temperature
    T2: outlet temperature
    eta: efficiency of the compressor
    gam: gamma value of the gases
    returns: pressure ratio
    """
    p1 = T2/T1
    p2 = gam/(gam-1)
    val = p1**(p2*eta)
    return val

def polypit(T1,T2,eta,gam): # For a turbine given temperatures
    """
    T1: inlet temperature
    T2: outlet temperature
    eta: efficiency of the turbine
    gam: gamma value of the gases
    returns: pressure ratio
    """
    p1 = T2/T1
    p2 = gam/(gam-1)
    val = p1**(p2/eta)
    return val

# Temperature ratio based on polytropic efficiency and pressures

def polytt(p1,p2,eta,gam): 
    """
    p1: inlet pressure
    p2: outlet pressure
    eta: efficiency of the turbine
    gam: gamma value of the gases
    returns: temperature ratio
    """
    p1 = p2/p1
    p2 = gam/(gam-1)
    val = p1**(p2*eta)
    return val

def polytc(pi,eta,gam):
    """
    pi: pressure ratio
    eta: efficiency of the compressor
    gam: gamma value of the gases
    returns: temperature ratio
    """
    p2 = (gam-1)/(gam)
    val = pi**(p2/eta)
    return val

# Temperatures after mixing

def mixingt(T1,T2,gam1,gam2,m1,m2):
    """
    T1: inlet temperature
    T2: inlet temperature
    gam1: gamma value of inlet gases at T1
    gam2: gamma value of inlet gases at T2
    m1: mass flow rate of inlet gases at T1
    m2: mass flow rate of inlet gases at T2
    returns: outlet temperature T3
    """
    cp1 = (gam1*R)/(gam1-1)
    cp2 = (gam2*R)/(gam2-1)
    m = m1/m2
    t1 = m*(cp1/cp2)
    t2 = t1*T1 - T2
    T3 = t2/(t1-1)
    return T3

# Effectiveness of the recuperator

def recup(e,Tnr,Tc_in,Th_in,mh,mc,gamh,gamc):
    """
    e: Effectiveness of the recuperator
    Tnr: known in the numerator - typically hot/cold inlet
    Tc_in: cold inlet
    Th_in: hot inlet
    mh: hot mass flow rate
    mc: cold mass flow rate
    gamh: gamma value of hot gases
    gamc: gamma value of cold gases
    returns: hot/cold outlet
    """
    m = mh/mc
    cph = (gamh*R)/(gamh-1)
    cpc = (gamc*R)/(gamc-1)
    cp = cph/cpc
    t1 = m*cp
    t2 = e/t1
    T = Tnr - t2*(Th_in-Tc_in)
    return T

# Heat exchange across a recuperator

def hxrecup(Tc_in,Th_in,e):
    """
    Tc_in: cold inlet
    Th_in: hot inlet
    returns: cold outlet
    """
    T = Tc_in + e*(Th_in-Tc_in)
    return T

# Effectiveness of the HX

def recir(e,Tnr,Tc_in,Th_in):
    """
    e: Effectiveness of the recirculator
    Tnr: known in the numerator - typically hot/cold inlet
    Tc_in: cold inlet
    Th_in: hot inlet
    returns: hot/cold outlet
    """
    T = Tnr - e*(Th_in-Tc_in)
    return T

# Local fuel-air ratio

def farat(Tin,Tout,gamh,gamc):
    """
    Tin: inlet temperature
    Tout: outlet temperature
    gamh: gamma value of the hot gases
    gamc: gamma value of cold gases
    returns: fuel-air ratio
    """
    gam = (gamh+gamc)/2
    cp = (gam*R)/(gam-1)
    t1 = q_f/cp
    f = (Tout-Tin)/(t1 - Tout)
    return f

# Turbocharger

def trunc(x,n):
    t = len('%.*x' % (n,x))
    return str(x)[:t]

def turboc(Tcin,Ttin,etat,etac,gamc,gamt):
    """
    Tcin: compressor inlet temperature
    Ttin: turbine inlet temperature
    etat: polytropic efficiency of the turbine
    etac: polytropic efficiency of the compressor
    gamc: gamma value of the compressor gases
    gamt: gamma value of the turbine gases
    returns: Turbocharger pressure ratio
    """
    p1 = (gamc-1)/gamc
    p2 = (gamt-1)/gamt
    T = Ttin/Tcin
    eta = etat*etac
    pi = 3
    x=1.001
    count = 0
    while x>1:
        # print "At T1 = "+str(Tcin)+" and T72 = "+str(Ttin)
        lhs = (x**(p1)) - 1
        # print "lhs = "+str(lhs)
        rhs = eta*T*(1-((1/x)**(p2)))
        # print "rhs = "+str(rhs)
        if (trunc(lhs,6) == trunc(rhs,6)):
            pi = x
            break
        elif count>1000000: # To avoid divergence
            break
        else:
            x=x+0.001
            count=count+1
    return pi

# The complete system

def system(T1,T32,T6,T7,T10,pi_cir,pi_cr,pi_hr,\
           pi_b,eff_cup,eff_cir,eta_tp,eta_cp,p1,size=1):
    """
    T1: inlet to LPC - same as ambient temperature
    T32: inlet to HPC
    T6: inlet to HPT
    T7: inlet to RC
    pi_cir: circulator pressure drop
    pi_hr: recuperator hot side pressure drop
    pi_cr: recuperator cold side pressure drop
    pi_b: burner pressure drop
    eff_cup: Effectiveness of recuperator
    eff_cir: Effectiveness of recirculator
    eta_tp: turbine polytropic efficiency
    eta_cp: compressor polytropic efficiency
    p1: ambient pressure - inlet to LPC
    returns: efficiency of the system and external cooling load per power
    """

    if size==1:
        print "For a large engine, \n"
    else:
        print "For a small engine, \n"

    # Finding pressures and temperatures
    pi_76 = polypit(T6,T7,eta_tp,gam_h)
    print "Pressure drop across the HPT is "+str(pi_76)
    pi_919 = 1
    pi_310 = 1
    pi_54 = pi_cr
    pi_65 = pi_b
    pi_97 = pi_hr
    pi_1091 = pi_cir
    pi_313 = pi_cir
    pi_3231 = pi_cir
    pi_324 = pi_3231*pi_313*pi_310*pi_1091*pi_919*pi_97*pi_76*pi_65*pi_54
    pi_432 = 1/pi_324
    print "Pressure ratio of the HPC is "+str(pi_432)
    tr_432 = polytc(pi_432,eta_cp,gam_m)
    T4 = tr_432*T32
    print "Temperature at the outlet of the HPC is "+str(T4)+" K"
    T5 = hxrecup(T4,T7,eff_cup)
    print "Temperature at the inlet of the burner is "+str(T5)+" K"
    f = farat(T5,T6,gam_h,gam_h)
    print "The local fuel air ratio is "+str(f)
    if size==1:
        r = (fg/(f*0.88))-1
        print "The recirculation ratio is "+str(r)
    else:
        r = (fg/f)-1
        print "The recirculation ratio is "+str(r)
    T9 = recup(eff_cup,T7,T4,T7,1+f,1,gam_h,gam_m)
    print "Temperature after recuperator is "+str(T9)+" K"
    T72 = T9
    pi_21 = turboc(T_amb,T72,eta_tp,eta_cp,gam_c,gam_h)
    print "The Turbocharger compression ratio is "+str(pi_21)
    tr_21 = polytc(pi_21,eta_cp,gam_c)
    T2 = tr_21*T1
    print "Temperature at the exit of LPT is "+str(T2)+" K"
    T3 = mixingt(T2,T10,gam_c,gam_m,1,r)
    print "Temperature after mixing is "+str(T3)+" K"
    p2 = pi_21*p1
    print "Pressure after LPC is "+str(p2)+" atm"
    p3 = p2
    p31 = pi_313*p3
    p32 = pi_3231*p31
    p4 = pi_432*p32
    p5 = pi_54*p4
    p6=pi_65*p5
    p7 = pi_76*p6
    p9 = pi_97*p7
    p91 = p9
    p72 = p91
    p8 = pi_21*p72
    print "Pressure after LPT is "+str(p8)+" atm"
    T91 = T9
    T72 = T91
    T31 = recir(eff_cir,T3,T_amb,T3)
    print "Temperature after the intercooler is "+str(T31)+" K"

    # Heat, work, power, etc
    q_gen = r*cp_h*(T91-T10)
    print "Heat by the generator is "+str(q_gen)
    q_eva = (1+r)*cp_m*(T31-T32)
    print "Heat by the evaporator is "+str(q_eva)
    q_ext = (cop*q_gen)-q_eva
    print "External heating load is "+str(q_ext)
    q_in = fg*q_f
    print "Input heat is "+str(q_in)
    if size==1:
        w = (0.88*(1+r))*(1+f)*cp_h*(T6-T7)-(1+r)*cp_m*(T4-T32)
    else:
        w = (1+r)*(1+f)*cp_h*(T6-T7)-(1+r)*cp_m*(T4-T32)
    eta = w/q_in
    print "The efficiency of the system is "+str(eta)
    q_cool = q_ext/w
    print "The external cooling load per power is "+str(q_cool)
    if size ==1:
        return (eta,q_cool,w)
    else:
        return (eta,q_cool,w,T4)

##########################################
# Analysis                               #
##########################################

# Design points

(eta_l,q_load_l,w_l) =\
        system(T_amb,T32d,T6_ld,T7d,T10d,pi_cir_ld,pi_crd,pi_hrd,\
               pi_b_ld,eff_cup_ld,eff_cir_ld,eta_tp_ld,eta_cp_ld,p_amb)

T4s_d = 0
(eta_s,q_load_s,w_s,T4s_d) = \
        system(T_amb,T32d,T6_sd,T7d,T10d,pi_cir_sd,pi_crd,pi_hrd,\
              pi_b_sd,eff_cup_sd,eff_cir_sd,eta_tp_sd,eta_cp_sd,p_amb,2)

########### Optimization ################

# Input parameters

T6l_opt = np.array([1600.0,1650.0,1700.0,1750.0])
T6s_opt = np.array([1200.0,1250.0,1300.0,1350.0])
T7_opt = np.array([1000.0,1050.0,1100.0,1150.0])
T32_opt = np.array([293.0,298.0,303.0,308.0])
T10_opt = np.array([368.0,373.0,378.0,383.0])
eta_opt_l = np.zeros(shape=(4,4,4,4))
q_opt_l = np.zeros(shape=(4,4,4,4))
w_opt_l = np.zeros(shape=(4,4,4,4))
eta_opt_s = np.zeros(shape=(4,4,4,4))
q_opt_s = np.zeros(shape=(4,4,4,4))
w_opt_s = np.zeros(shape=(4,4,4,4))

# Large engine

# for i in range(0,4):
    # for j in range(0,4):
        # for k in range(0,4):
            # for l in range(0,4):
                # print "\nFor T6 = "+str(T6l_opt[l])+", T7 = "+str(T7_opt[k])+\
                        # ", T10 = "+str(T10_opt[j])+" T32 = "+str(T32_opt[i])+","
                # (eta_opt_l[i,j,k,l],q_opt_l[i,j,k,l],w_opt_l[i,j,k,l]) = \
                    # system(T_amb,T32_opt[i],T6l_opt[l],T7_opt[k],T10_opt[j],\
                         # pi_cir_ld,pi_crd,pi_hrd,pi_b_ld,eff_cup_ld,eff_cir_ld,\
                           # eta_tp_ld,eta_cp_ld,p_amb)
                # print "\n"

eta_max_l = 0.5
imax_l = 0
jmax_l = 0
kmax_l = 0
lmax_l = 0

# for i in range(0,4):
    # for j in range(0,4):
        # for k in range(0,4):
            # for l in range(0,4):
                # if eta_opt_l[i,j,k,l]>eta_max_l:
                    # eta_max_l = eta_opt_l[i,j,k,l]
                    # imax_l = i
                    # jmax_l = j
                    # kmax_l = k
                    # lmax_l = l

# print "\nThe maximum efficiency of "+str(eta_max_l)+\
        # " for the large engine when optimizing occurs at T32 = "\
        # +str(T32_opt[imax_l])+", T10 = "+str(T10_opt[jmax_l])+", T7 = "+\
        # str(T7_opt[kmax_l])+" and T6 = "+str(T6l_opt[lmax_l])+\
        # "\nwith a cooling load of "+str(q_opt_l[imax_l,jmax_l,kmax_l,lmax_l])+\
        # " and specific work of "+str(w_opt_l[imax_l,jmax_l,kmax_l,lmax_l])+\
        # "\nwhile the efficiency at Design point is "+str(eta_l)+"\n"

# Small engine

T4_opt_s = np.zeros(shape=(4,4,4,4))

# for i in range(0,4):
    # for j in range(0,4):
        # for k in range(0,4):
            # for l in range(0,4):
                # print "\nFor T6 = "+str(T6l_opt[l])+", T7 = "+str(T7_opt[k])+\
                        # ", T10 = "+str(T10_opt[j])+" T32 = "+str(T32_opt[i])+","
                # (eta_opt_s[i,j,k,l],q_opt_s[i,j,k,l],w_opt_s[i,j,k,l],\
                 # T4_opt_s[i,j,k,l]) = system(T_amb,T32_opt[i],T6s_opt[l],T7_opt[k],\
                    # T10_opt[j],pi_cir_ld,pi_crd,pi_hrd,pi_b_ld,eff_cup_ld,eff_cir_ld,\
                        # eta_tp_ld,eta_cp_ld,p_amb,2)
                # print "\n"

eta_max_s = 0.5
imax_s = 0
jmax_s = 0
kmax_s = 0
lmax_s = 0
                
# for i in range(0,4):
    # for j in range(0,4):
        # for k in range(0,4):
            # for l in range(0,4):
                # if eta_opt_s[i,j,k,l]>eta_max_s:
                    # eta_max_s = eta_opt_s[i,j,k,l]
                    # imax_s = i
                    # jmax_s = j
                    # kmax_s = k
                    # lmax_s = l

# print "\nThe maximum efficiency of "+str(eta_max_s)+\
        # " for the small engine when optimizing occurs at T32 = "\
        # +str(T32_opt[imax_s])+", T10 = "+str(T10_opt[jmax_s])+", T7 = "+\
        # str(T7_opt[kmax_s])+", T6 = "+str(T6l_opt[lmax_s])+" and T4 = "+\
        # str(T4_opt_s[imax_s,jmax_s,kmax_s,lmax_s])+\
        # "\nwith a cooling load of "+str(q_opt_s[imax_s,jmax_s,kmax_s,lmax_s])+\
        # " and specific work of "+str(w_opt_s[imax_s,jmax_s,kmax_s,lmax_s])+\
        # "\nwhile the efficiency at Design point is "+str(eta_s)+"\n"

##############################################################################
# The End                                                                    #
##############################################################################

