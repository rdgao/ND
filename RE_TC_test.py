from brian2 import *
import scipy as sp
import pylab as plt

plt.close('all')
seed(0)
# Parameters
Cm = 1.*ufarad


###########################RE#################################
E_NaRE = 50.*mV
g_NaRE = 200.*msiemens
E_KRE = -100.*mV
g_KRE = 20.*msiemens
E_CaRE = 120.*mV
g_CaRE = 3.*msiemens
E_LRE = -90.*mV
g_LRE = 0.05*msiemens

model_RE = Equations('''
    
Vt_Na = v + 55.*mV : volt
Vt_K = v + 55.*mV  : volt

am = ( (0.32/mV)*(13.*mV - Vt_Na) / (exp((13.*mV - Vt_Na)/(4.*mV))-1.) ) /ms : Hz
bm = ( (0.28/mV)*(Vt_Na - 40.*mV) / (exp((Vt_Na - 40.*mV)/(5.*mV))-1.) )/ms : Hz
ah = 0.128 * exp((17.*mV-Vt_Na)/(18.*mV))/ms : Hz
bh = 4.0/(1. + exp((40.*mV-Vt_Na) / (5.*mV)))/ms : Hz 

an = (.032/mV*(15.*mV-Vt_K))/(exp((15.*mV-Vt_K)/(5.*mV))-1.)/ms : Hz
bn = ((0.5)*exp((10.*mV - Vt_K)/(40.*mV)))/ms : Hz


mT_inf = 1. / (1. + exp(-(v + 52.*mV)/(7.4*mV))) : 1
hT_inf = 1. / (1. + exp( (v + 80.*mV)/(5.*mV)))  : 1

tau_mT = (0.44 + 0.15/(exp((v+27.*mV)/(10.*mV))+exp(-(v+102.*mV)/(15.*mV)))) *ms : second
tau_hT = (22.7 + 0.27/(exp((v+48.*mV)/(4.*mV))+exp(-(v+407.*mV)/(50.*mV)))) *ms : second        

I_Na = g_NaRE * m**3 * h * (v - E_NaRE) : amp
I_K = g_KRE * n**4 * (v - E_KRE) : amp    
I_T = g_CaRE * mT**2 * hT * (v - E_CaRE) : amp
I_L = g_LRE * (v - E_LRE) : amp
I_app : amp

dm/dt = (am)*(1-m) - (bm)*m : 1
dh/dt = (ah)*(1-h) - (bh)*h : 1
dn/dt = (an)*(1-n) - (bn)*n : 1
dmT/dt = -(mT - mT_inf)/tau_mT : 1
dhT/dt = -(hT - hT_inf)/tau_hT : 1

dv/dt = (I_app - I_Na - I_K - I_T - I_L)/Cm : volt

''')

###########################TC#################################
E_NaTC = 50.*mV
g_NaTC = 90.*msiemens
E_KTC = -100.*mV
g_KTC = 10.*msiemens
E_CaTC = 120.*mV
g_CaTC = 3.*msiemens
E_hTC = -40.*mV
g_hTC = 0.025*msiemens


model_TC = Equations('''
    
Vt_Na = v + 35.*mV : volt
Vt_K = v + 25.*mV  : volt
Vt_T = v + 2.*mV : volt

### sodium and potassium ###
am = ( (0.32/mV)*(13.*mV - Vt_Na) / (exp((13.*mV - Vt_Na)/(4.*mV))-1.) ) /ms : Hz
bm = ( (0.28/mV)*(Vt_Na - 40.*mV) / (exp((Vt_Na - 40.*mV)/(5.*mV))-1.) )/ms : Hz
ah = 0.128 * exp((17.*mV-Vt_Na)/(18.*mV))/ms : Hz
bh = 4.0/(1. + exp((40.*mV-Vt_Na) / (5.*mV)))/ms : Hz 
an = (.032/mV*(15.*mV-Vt_K))/(exp((15.*mV-Vt_K)/(5.*mV))-1.)/ms : Hz
bn = ((0.5)*exp((10.*mV - Vt_K)/(40.*mV)))/ms : Hz

### T-current ####
mT_inf = 1. / (1. + exp(-(Vt_T + 57.*mV)/(6.2*mV))) : 1
tau_hT = ((30.8 + 
            (211.4 + exp((Vt_T+113.2*mV)/(5.*mV)))
            /(1.+exp((Vt_T+84.*mV)/(3.2*mV)))
            ) /3.73) *ms :second
hT_inf = 1. / (1. + exp( (Vt_T + 81.*mV)/(4.*mV)))  : 1


### h-current ####
hTC_inf = 1./(1.+exp((v+75.*mV)/(5.5*mV))) : 1
tau_s = (20. + 1000./(
            (exp((v+71.5*mV)/(14.2*mV)))
            + (exp(-(v+89.*mV)/(11.6*mV)))
            ))*ms : second

a_hTC = hTC_inf/tau_s :Hz
b_hTC = (1.-hTC_inf)/tau_s :Hz


I_Na = g_NaTC * m**3 * h * (v - E_NaTC) : amp
I_K = g_KTC * n**4 * (v - E_KTC) : amp    
I_T = g_CaTC * mT_inf * hT * (v - E_CaTC) : amp
I_h = g_hTC*(o_1 + 2.*(1.-c_1-o_1))*(v-E_hTC) :amp
I_L = 0.01*msiemens * (v - (-70.*mV)) + 0.0172*msiemens *(v - (-100.*mV)) :amp
I_app : amp

dm/dt = (am)*(1-m) - (bm)*m : 1
dh/dt = (ah)*(1-h) - (bh)*h : 1
dn/dt = (an)*(1-n) - (bn)*n : 1
dhT/dt = -(hT - hT_inf)/tau_hT : 1

dc_1/dt = b_hTC*o_1-a_hTC*c_1 : 1
dp_0/dt = (0.0004*(1.-p_0)-0.0004*(C_CaInt/0.002)**4 * p_0)/ms : 1
do_1/dt = (0.001*(1.-c_1-o_1) - 0.001*((1.- p_0)/0.01)*o_1)/ms : 1

It1 = ( ((-10.*(I_T/amp)/(2.*96489.))>0)*(-10.*(I_T/amp)/(2.*96489.)) )/ms :Hz
It3 = ((0.00024 - C_CaInt)/5.0)/ms :Hz
dC_CaInt/dt = It1 + It3 : 1

dv/dt = (I_app - I_Na - I_K - I_T - I_h - I_L)/Cm : volt

''')





dt_sim = 0.1*ms
duration = 100.*ms

num_neurons = 1
N = NeuronGroup(num_neurons, model_TC,
                    threshold='not_refractory and (v > 0*mV)',
                    refractory='v > 0*mV')
N.v = -70.*mV                    
m = StateMonitor(N, ['v', 'I_Na', 'I_K', 'I_T', 'I_h', 'I_L'], record = True, dt = dt_sim)
m_s = SpikeMonitor(N)

N.I_app = 20.*uA# 5.*uA * N.i / num_neurons

run(duration)

figure()
plot(N.I_app/nA, m_s.count / duration)
xlabel('I (nA)')
ylabel('Firing rate (sp/s)')
show()

figure()
plot(m.t, m.v[-1])
show()