from brian2 import *
import scipy as sp
import pylab as plt


#TODO as of 11/24
#1. Correct number of neurons (80PY; 16LTS,FS; 6 RE,TC)=
#2. Add tonic current to PY (1.8), FS (.55), and LTS (1). +normal(0,.1 var)
#3. Simulate EEG
#4. Make sure the simulation is behaving properly

plt.close('all')
seed(0)
# Parameters
Cm = 1.*ufarad

######################---pyramidal---#######################
#leak
g_L = 0.1*msiemens
E_L = -67*mV
#sodium
g_Na = 100*msiemens
E_Na = 50*mV
#potassium
g_K = 80*msiemens
E_K = -100*mV

model_PY = Equations('''

am = (0.32/mV) * (v+54*mV) / (1- exp(-(v+54*mV) / (4*mV)))/ms : Hz
bm = (0.28/mV) * (v+27*mV) / (exp((v+27*mV) / (5*mV))-1)/ms : Hz
ah = 0.128 * exp(-(v+50.*mV)/(18.*mV))/ms : Hz
bh = 4.0/(1. + exp(-(v+27.*mV) / (5.*mV)))/ms : Hz
an = (.032/mV*(v+52.*mV))/(1- exp(-(v+52.*mV)/(5.*mV)))/ms : Hz
bn = 0.5*exp(-(v+57.*mV)/(40.*mV))/ms : Hz

I_Na = g_Na * (m**3) * h * (v-E_Na) :amp
I_K = g_K * (n**4) * (v-E_K) :amp
I_L = g_L * (v-E_L) :amp

dn/dt = (an)*(1-n) - (bn)*n : 1
dm/dt = (am)*(1-m) - (bm)*m : 1
dh/dt = (ah)*(1-h) - (bh)*h : 1

dv/dt = (-I_L - I_Na - I_K + Iapp + IFS + ILTS + ITC + IPY)/Cm : volt
            
Iapp : amp
IPY : amp
IFS : amp
ILTS : amp
ITC : amp
''')

#Add IPY for the FS cells

#########################---LTS---################################
#M-current
g_m2 = 2*msiemens
E_m2 = -100*mV
Qs = 3.209
model_LTS = Equations('''

am = (0.32/mV) * (v+54*mV) / (1- exp(-(v+54*mV) / (4*mV)))/ms : Hz
bm = (0.28/mV) * (v+27*mV) / (exp((v+27*mV) / (5*mV))-1)/ms : Hz
ah = 0.128 * exp(-(v+50.*mV)/(18.*mV))/ms : Hz
bh = 4.0/(1. + exp(-(v+27.*mV) / (5.*mV)))/ms : Hz
an = (.032/mV*(v+52.*mV))/(1- exp(-(v+52.*mV)/(5.*mV)))/ms : Hz
bn = 0.5*exp(-(v+57.*mV)/(40.*mV))/ms : Hz
am2 = 0.0001/mV * Qs * (v+30.*mV) / (1 - exp(-(v+30.*mV)/(9.*mV)))/ms : Hz
bm2 =-0.0001/mV * Qs * (v+30.*mV) / (1 - exp((v+30.*mV)/(9.*mV)))/ms : Hz

I_Na = g_Na*(m**3)*h*(v-E_Na) :amp
I_K = g_K*(n**4)*(v-E_K) :amp
I_Ca = g_m2*m2*(v-E_m2) :amp
I_L = g_L*(v-E_L) :amp

dn/dt = (an)*(1-n) - (bn)*n : 1
dm/dt = (am)*(1-m) - (bm)*m : 1
dh/dt = (ah)*(1-h) - (bh)*h : 1
dm2/dt = (am2)*(1-m2) - (bm2)*m2 : 1

dv/dt = (-I_L - I_Na - I_K - I_Ca + Iapp + IPY + ILTS + ITC)/Cm : volt

Iapp : amp
IPY : amp
ILTS : amp
ITC : amp
''')

#######################----FS----###############################
#FS cells use same model as PY cells, with different applied current

########################----TC----##############################
#sodium
g_NaTC = 90.*msiemens
E_NaTC = 50.*mV
#potassium
g_KTC = 10.*msiemens
E_KTC = -100.*mV
#T-current
g_CaTC = 2.*msiemens
E_CaTC = 120.*mV
#h-current
g_hTC = 0.25*msiemens
E_hTC = -40.*mV

model_TC = Equations('''

Vtm = v + 35.*mV : volt
Vtn = v + 25.*mV : volt
Vth2= v + 2.*mV : volt

an = (.032/mV*(15.*mV-Vtn))/(exp((15.*mV-Vtn)/(5.*mV))-1.)/ms : Hz
bn = (((0.5/mV) * (10.*mV - Vtn)) / (exp((10.*mV-Vtn) / (40.*mV))))/ms : Hz

am = (0.32/mV) * (13.*mV - Vtm) / (exp((13.*mV-Vtm) / (4.*mV)) -1.)/ms : Hz
bm = (0.28/mV) * (Vtm - 40.*mV) / (exp((Vtm - 40.*mV) / (5.*mV))-1.)/ms : Hz

ah = 0.128 * ((17.*mV-Vtm)/(18.*mV))/ms : Hz
#ah = 0.128 * exp((17.*mV-Vtm)/(18.*mV))/ms : Hz
bh = 4.0/(1. + exp((40.*mV-Vtm) / (5.*mV)))/ms : Hz

h2_inf = 1. / (1. + exp((Vth2 + 81.*mV)/(4.*mV))) : 1
m2_inf =  1. / (1. + exp(-(Vth2 + 57.*mV)/(6.2*mV))) : 1

tau_h2 = ((30.8 + 
            (211.4 + exp((Vth2+113.2*mV)/(5.*mV)))
            /(1.+exp((Vth2+84.*mV)/(3.2*mV)))
            ) /3.73) *ms :second

a_hTC = hTC_inf/tau_s :Hz
b_hTC = (1.-hTC_inf)/tau_s :Hz

hTC_inf = 1./(1.+exp((v+75.*mV)/(5.5*mV))) : 1
tau_s = (20. + 1000./(
            (1.+exp((v+71.5*mV)/(14.2*mV)))
            + (1. +exp(-(v+89.*mV)/(11.6*mV)))
            ))*ms : second

dC_CaInt/dt = (-10.*(I_T/amp)/(2.*96489.) + (0.00024 - C_CaInt)/5.)*Hz : 1

I_L = 0.01*msiemens * (v - (-70.*mV)) + 0.0172*msiemens *(v - (-100.*mV))  :amp
I_Na = g_NaTC * (m**3) * h *(v-E_NaTC) :amp
I_K = g_KTC * (n**4) * (v-E_KTC) :amp
I_h = g_hTC*(o_1 + 2.*(1.-c_1-o_1))*(v-E_hTC) :amp
I_T = g_CaTC*(m2_inf**2)*h2*(v-E_CaTC) :amp

dv/dt = (-I_L - I_Na - I_K - I_T - I_h + IPY + IRE)/Cm : volt
        
dn/dt = (an)*(1.-n) - (bn)*n : 1
dm/dt = (am)*(1.-m) - (bm)*m : 1
dh/dt = (ah)*(1.-h) - (bh)*h : 1
dh2/dt= (h2_inf - h2) / tau_h2 : 1
do_1/dt = (0.001*(1.-c_1-o_1)-0.001*((1.- p_0)/0.01))*Hz : 1
dp_0/dt = (0.0004*(1.-p_0)-0.0004*(C_CaInt/0.002)**4)*Hz : 1
dc_1/dt = b_hTC*o_1-a_hTC*c_1 : 1

IPY : amp
IRE : amp

''')


########################----RE----##############################
#sodium
g_NaRE = 200.*msiemens
E_NaRE = 50.*mV
#potassium
g_KRE = 20.*msiemens
E_KRE = -100.*mV
#T-current
g_CaRE = 3.*msiemens
E_CaRE = 120.*mV
#leak
g_LRE = 0.05*msiemens
E_LRE = -90.*mV

g_hRE = 0.25*msiemens

model_RE = Equations('''

Vtm = v + 55.*mV : volt
Vtn = v + 55.*mV : volt
Vth2= v + 2.*mV : volt

an = (.032/mV*(15.*mV-Vtn))/(exp((15.*mV-Vtn)/(5.*mV))-1.)/ms : Hz
bn = (((0.5/mV) * (10.*mV - Vtn)) / (exp((10.*mV-Vtn) / (40.*mV))))/ms : Hz

am = (0.32/mV) * (13.*mV - Vtm) / (exp((13.*mV-Vtm) / (4.*mV)) -1.)/ms : Hz
bm = (0.28/mV) * (Vtm - 40.*mV) / (exp((Vtm - 40.*mV) / (5.*mV))-1.)/ms : Hz
ah = 0.128 * ((17.*mV-Vtm)/(18.*mV))/ms : Hz
#ah = 0.128 * exp((17.*mV-Vtm)/(18.*mV))/ms : Hz

bh = 4.0/(1. + exp((40.*mV-Vtm) / (5.*mV)))/ms : Hz

h2_inf = 1. / (1. + exp((v + 80.*mV)/(5.*mV))) : 1
m2_inf =  1. / (1. + exp(-(v + 52.*mV)/(7.4*mV))) : 1

tau_h2 = (22.7 + 0.27/(exp((v+48.*mV)/(4.*mV))+exp(-(v+407.*mV)/(50.*mV))))*ms : second
tau_m2 = (0.44 + 0.15/(exp((v+27.*mV)/(10.*mV))+exp(-(v+102.*mV)/(15.*mV))))*ms : second

I_Na = g_NaRE*(m**3)*h*(v-E_NaRE) :amp
I_K = g_KRE*(n**4)*(v-E_KRE) :amp
I_T = g_CaRE * (m2**2) * h2 *(v-E_CaRE) :amp
I_L = g_LRE*(v-E_LRE) :amp

dn/dt = (an)*(1-n) - (bn)*n : 1
dm/dt = (am)*(1-m) - (bm)*m : 1
dh/dt = (ah)*(1-h) - (bh)*h : 1
dh2/dt= (h2_inf - h2) / tau_h2 : 1
dm2/dt= (m2_inf - m2) / tau_m2 : 1

dv/dt = (-I_L - I_Na - I_K - I_T + IPY + ITC + IRE)/Cm : volt
   
IPY : amp
ITC : amp
IRE : amp

''')

###########################
#CREATE neuron groups
inp  = PoissonGroup(1, 100*Hz)
PY = NeuronGroup(2, model_PY,threshold='not_refractory and (v > 0*mV)',
                    refractory='v > 0*mV')
LTS =NeuronGroup(2, model_LTS,threshold='not_refractory and (v > 0*mV)',
                    refractory='v > 0*mV')
FS = NeuronGroup(2, model_PY,threshold='not_refractory and (v > 0*mV)',
                    refractory='v > 0*mV')
TC = NeuronGroup(2, model_TC,threshold='not_refractory and (v > 0*mV)',
                    refractory='v > 0*mV')
RE = NeuronGroup(2, model_RE,threshold='not_refractory and (v > -40.*mV)',
                    refractory='v > -40.*mV')
PY.v = E_L
FS.v = E_L
LTS.v = E_L
TC.v = -85.*mV
RE.v = E_LRE  

#PY.Iapp = 1.8
#LTS.Iapp= 1
#FS.Iapp = .55                

##########################
#CREATE synapses
#15 synapses
# Poisson to PY
# PY to: FS, LTS, TC, RE
# FS to: PY, FS
# LTSto: PY, LTS
# RE to: TC, RE
# TC to: FS, LTS, PY, RE

E_AMPA = 0*mV
E_GABA = -80*mV
tau_GABA = .005 #5ms
tau_AMPA = .002

g_AMPA_PY = .1*siemens
g_AMPA_LTS = .5*siemens
g_AMPA_FS = 2*siemens
g_AMPA_RE = .1*siemens
g_AMPA_TC = .1*siemens
g_GABA_PY = .64*siemens
g_GABA_LTS = .15*siemens
g_GABA_FS = 1*siemens
g_GABA_RE = .06*siemens
g_GABA_TC = .06*siemens

c_inpPY       = Synapses(inp, PY,model='w : volt', pre='v += w')
c_inpPY.connect(True) #all-to-all connectivity
c_inpPY.w     = '50. * mV'
c_inpPY.delay = 2*ms 

c_PYLTS       = Synapses(PY, LTS,
                         model='''IPY_post = -g_AMPA_LTS*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_PYLTS.connect(True)
c_PYLTS.delay = 5*ms 

c_LTSPY       = Synapses(LTS, PY,
                         model='''ILTS_post = -g_GABA_PY*x*(v_post-E_GABA) : amp (summed)
                         dx/dt = 2 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_GABA : 1
                         ''')
c_LTSPY.connect(True)
c_LTSPY.delay = 5*ms 

c_PYFS       = Synapses(PY, FS,
                         model='''IPY_post = -g_AMPA_FS*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_PYFS.connect(True)
c_PYFS.delay = 5*ms 

c_FSPY       = Synapses(FS, PY,
                         model='''IFS_post = -g_GABA_PY*x*(v_post-E_GABA) : amp (summed)
                         dx/dt = 2 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_GABA : 1
                         ''')
c_FSPY.connect(True)
c_FSPY.delay = 5*ms 

c_PYTC       = Synapses(PY, TC,
                         model='''IPY_post = -g_AMPA_TC*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_PYTC.connect(True)
c_PYTC.delay = 5*ms 

c_TCPY       = Synapses(TC, PY,
                         model='''ITC_post = -g_AMPA_PY*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_TCPY.connect(True)
c_TCPY.delay = 5*ms 

c_TCFS       = Synapses(TC, FS,
                         model='''ITC_post = -g_AMPA_FS*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_TCFS.connect(True)
c_TCFS.delay = 5*ms 

c_TCLTS       = Synapses(TC, LTS,
                         model='''ITC_post = -g_AMPA_LTS*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_TCLTS.connect(True)
c_TCLTS.delay = 5*ms 

c_PYRE       = Synapses(PY, RE,
                         model='''IPY_post = -g_AMPA_RE*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_PYRE.connect(True)
c_PYRE.delay = 5*ms 

c_TCRE       = Synapses(TC, RE,
                         model='''ITC_post = -g_AMPA_RE*x*(v_post-E_AMPA) : amp (summed)
                         dx/dt = 5 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_AMPA : 1
                         ''')
c_TCRE.connect(True)
c_TCRE.delay = 5*ms 

c_RETC       = Synapses(RE, TC,
                         model='''IRE_post = -g_GABA_TC*x*(v_post-E_GABA) : amp (summed)
                         dx/dt = 2 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_GABA : 1
                         ''')
c_RETC.connect(True)
c_RETC.delay = 5*ms 

c_RE       = Synapses(RE, model='''IRE_post = -g_GABA_RE*x*(v_post-E_GABA) : amp (summed)
                         dx/dt = 2 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_GABA : 1
                         ''')
c_RE.connect(True)
c_RE.delay = 1*ms 

c_FS       = Synapses(FS, model='''IFS_post = -g_GABA_FS*x*(v_post-E_GABA) : amp (summed)
                         dx/dt = 2 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_GABA : 1
                         ''')
c_FS.connect(True)
c_FS.delay = 1*ms 

c_LTS       = Synapses(LTS, model='''ILTS_post = -g_GABA_LTS*x*(v_post-E_GABA) : amp (summed)
                         dx/dt = 2 * (1 + tanh(v_pre/(4*mV))) * (1-x) - x/tau_GABA : 1
                         ''')
c_LTS.connect(True)
c_LTS.delay = 1*ms 


#####################
#Define monitors for the simulation
dt_sim = 0.1*ms

m_PY      = StateMonitor(PY, 'v', record = True, dt = dt_sim)
m_LTS     = StateMonitor(LTS,'v', record = True, dt = dt_sim)
m_FS      = StateMonitor(FS, 'v', record = True, dt = dt_sim)
m_TC      = StateMonitor(TC, 'v', record = True, dt = dt_sim)
m_RE      = StateMonitor(RE, 'v', record = True, dt = dt_sim)
ms_PY   = SpikeMonitor(PY)
ms_inp = SpikeMonitor(inp)

s_PYLTS      = StateMonitor(c_PYLTS,'x', record = True)
s_LTSPY      = StateMonitor(c_LTSPY,'x', record = True)
I_LTSPY      = StateMonitor(PY,'ILTS', record = True, dt = dt_sim)
I_PYLTS      = StateMonitor(LTS,'IPY', record = True, dt = dt_sim)
I_FSPY       = StateMonitor(PY,'IFS', record = True, dt = dt_sim)
I_PYFS       = StateMonitor(FS,'IPY', record = True, dt = dt_sim)

#########################
# Run the simulation
monitors = [m_PY,m_LTS,m_FS,ms_inp,#m_TC,m_RE,
            s_PYLTS,s_LTSPY,I_LTSPY,I_PYLTS,I_FSPY,I_PYFS]
connects = [c_inpPY, c_PYLTS, c_LTSPY, c_PYFS, c_FSPY]#,
            #c_PYTC, c_TCPY, c_TCFS, c_TCLTS,
            #c_PYRE, c_TCRE, c_RETC, c_RE, c_FS, c_LTS]
neurons = [inp, PY, LTS, FS]#, TC, RE]
net = Network(neurons,connects,monitors)
net.run(100*ms)


#########################
# Visualize results
#%%
plt.figure()
subplot(3,1,1)
plt.plot(m_PY.t/ms, m_PY.v[0]/mV, label='PY')
plt.plot(m_LTS.t/ms, m_LTS.v[0]/mV, label='LTS')
plt.plot(m_FS.t/ms, m_FS.v[0]/mV, label='FS')
plt.plot(m_TC.t/ms, m_TC.v[0]/mV, label='TC')
plt.plot(m_RE.t/ms, m_RE.v[0]/mV, label='RE')
plt.plot(ms_inp.t/ms, ms_inp.i, 'ro')
plt.legend()
plt.show()
#%%
subplot(3,1,2)
plt.plot(s_PYLTS.t/ms,s_PYLTS.x[0],label='PYLTS')
plt.plot(s_LTSPY.t/ms,s_LTSPY.x[0],label='LTSPY')
plt.legend()
plt.show()
subplot(3,1,3)
plt.plot(I_LTSPY.t/ms,I_LTSPY.ILTS[0],label='I_LTS-->PY')
plt.plot(I_PYLTS.t/ms,I_PYLTS.IPY[0],label='I_PY-->LTS')
plt.plot(I_FSPY.t/ms,I_FSPY.IFS[0],label='I_FS-->PY')
plt.plot(I_PYFS.t/ms,I_PYFS.IPY[0],label='I_PY-->FS')
plt.legend()
plt.show()