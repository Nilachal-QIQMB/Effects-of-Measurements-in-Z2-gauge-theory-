#!/usr/bin/env python
# coding: utf-8

# In[2]:


pip install qutip


# In[3]:


from qutip import *


# In[6]:


import numpy as np
import matplotlib.pyplot as plt


# In[7]:


def pauliop(mat,nmax,i):
  pauliop_mat= [qeye(2)]*nmax
  spinop= mat
  pauliop_mat[i]=spinop
  return (tensor(pauliop_mat))


# In[8]:


sigmaz()


# In[9]:


pauliop(sigmaz(),4,0)


# In[10]:


pauliop(qeye(2),4,0)


# In[11]:


def numop(nmax,i):
    idmat= pauliop(qeye(2),nmax,i)
    spinz= pauliop(sigmaz(),nmax,i)
    numop= 0.5*(idmat-spinz)
    return numop


# In[12]:


def total_numop(nmax):
    num_op=0
    for i in range(nmax):
     if(i%2==0):
      num_op= num_op+numop(nmax,i)
    return num_op


# In[13]:


def Z2_Ham_nh_new1(nsites,X,f,mu,gamma):
 vertex= nsites 
 links= nsites-1
 ndim= vertex+links
 b=(ndim/2)
 ## Constructing the Hamiltonian
 H_J=0
 H_f=0
 H_m=0
 for i in range(0,ndim-2,2):
  t1= pauliop(sigmap(),ndim,i)*0.5*pauliop(sigmax(),ndim,i+1)*pauliop(sigmam(),ndim,i+2)
  H_J+= (X)*(t1+t1.dag())
  #if 0<= i < b-2 :
   #H_J+= (X-(1j)*gamma)*(t1+t1.dag())
  #else:
    #H_J+= (X)*(t1+t1.dag())
 for j in range(1,ndim-1,2):
  H_f+= (f-(1j*gamma))*0.5*pauliop(sigmaz(),ndim,j)
 for k in range(0,ndim,2):
  vertex_index = (k) // 2
  H_m += mu * ((-1)**vertex_index)*0.5*pauliop(sigmaz(),ndim,k)

 Full_Ham= H_J+ H_f + H_m
 return(Full_Ham)


# In[14]:


from qutip import basis, tensor

def strong_coupling_vacuum(nsites):
    """
    Gauge-invariant strong coupling vacuum for Z2 LGT
    Ordering: v0, l0, v1, l1, ..., v_{N-1}
    """
    state_list = []

    for i in range(0,(2*nsites-1)):
        if i % 2 == 0:
           # vertex (matter)
            k = (i // 2)
            if k % 2 == 0:
                state_list.append(basis(2, 1))  # |down>
            else:
                state_list.append(basis(2, 0))  # |up>
        else:
            # link (gauge)
            state_list.append(basis(2, 1))      # |down>

    return tensor(state_list)


# In[15]:


def Electric_term(f,dim):
 HE=0
 for j in range(1, dim-1, 2):
   Electric_term = f* 0.5*pauliop(sigmaz(),dim, j)
   HE+= Electric_term
 return HE


# In[16]:


nsites= 6
dim=(2*nsites-1)
H_herm1= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0)
psi0 = strong_coupling_vacuum(nsites)
tvals=  np.arange(0,10.1,0.1)
f=1.0
eops= [Electric_term(f,dim)]
result= sesolve(H_herm1,psi0,tvals,eops)
print(result.expect[0])


# In[17]:


eiglist1= H_herm1.eigenenergies()
print(eiglist1)
for eig  in eiglist1:
 print(eig)


# In[18]:


result_without_state = sesolve(H_herm1,psi0,tvals)
states= result_without_state.states
psi0 = strong_coupling_vacuum(nsites)
Pt_full_evolution=[]
for state in states:
  overlap= psi0.dag()*state
  p_t= np.linalg.norm(overlap)
  print(p_t**2)
  Pt_full_evolution.append(p_t**2)


# In[11]:


## plotting
check1= np.loadtxt("data_check_1storder.dat")
check2= np.loadtxt("datacheck_2ndorder.dat")
check3= np.loadtxt("data_check_TDVP.dat")

## x-list and y_list
x_check1= check1[:,0]
y_check1= check1[:,1]

x_check2= check2[:,0]
y_check2= check2[:,1]

x_check3= check3[:,0]
y_check3= check3[:,1] 
## Plotting 
plt.figure(figsize=(12,6))
plt.plot(x_check1,y_check1,marker=".",label=r"$TEBD$ $1st$ $order$")
plt.plot(x_check2,y_check2,marker=".",label=r"$TEBD$ $2nd$ $order$")
plt.plot(x_check3,y_check3,marker=".",label=r"$TDVP$")
plt.plot(tvals,result.expect[0],marker=".",label="$Exact$")
plt.xlabel("$t$",fontsize=15)
plt.ylabel(r"$\langle E \rangle (t)$",fontsize=15)
plt.title("$\gamma=0,Hermitian$")
plt.legend()
plt.show()


# ## Gauss Law check at each site ##

# In[19]:


def gauss_law(ndim,i):
    index=i//2
    argument= (0.5*(pauliop(sigmaz(),ndim,i)+1)- 0.5*(1-((-1)**index)))
    phase= (-1j*np.pi*argument).expm()
    if (i==0):
      Gi= -1*pauliop(sigmaz(),ndim,i+1)* phase
    elif(i==ndim-1):
      Gi= pauliop(sigmaz(),ndim,i-1)*-1* phase
    else:
     Gi= pauliop(sigmaz(),ndim,i-1)*pauliop(sigmaz(),ndim,i+1)* phase
    return Gi


# In[20]:


import itertools

def all_possible_states(ndim):
 up = basis(2, 0)
 dn = basis(2, 1)

 basis_states = []

 for config in itertools.product([up, dn], repeat=ndim):
    psi = tensor(config)
    basis_states.append(psi)
 return basis_states


# In[21]:


all_possible_states(3)[0]


# In[22]:


nsites=6
links= nsites-1
ndim= nsites+links

possible_state_list= all_possible_states(ndim)
print(len(possible_state_list))


# In[23]:


psi0= strong_coupling_vacuum(nsites)
psi0.norm()


# In[24]:


for i in range(ndim):
 if(i< ndim-1 and i%2==0):
   Gi_exp= psi0.dag()*gauss_law(ndim,i)*psi0
   print(i,Gi_exp)


# In[25]:


for i in range(0,ndim,2):
#if(i%2==0):
 Gi_exp= psi0.dag()*gauss_law(ndim,i)*psi0
 print(i,Gi_exp)


# In[26]:


product=1
for i in range(0,ndim,2):
    term= gauss_law(ndim,i)
    product=product*term
print(product)


# ## Gauss law invariant subspace formation ##

# ## Initial calculation with N=4 ##

# In[27]:


#physical_state1=[]
nsites=4
links= nsites-1
ndim= nsites+links
states_list= all_possible_states(ndim)
print(len(states_list))
gauss_ops = {i: gauss_law(ndim,i)for i in range(0,ndim-2,2)}
val_master_list=[]
for state in states_list:
    val_list= np.zeros(len(gauss_ops),dtype=complex)
    for i, Gi in gauss_ops.items():
        val = (state.dag()*Gi*state)
        index=(i//2)
        val_list[index]=val
    val_master_list.append((val_list))


# In[28]:


my_list = [1, 1, 1, 1]

is_all_1= all(x == 1 for x in my_list)
print(is_all_1)
if is_all_1==True:
    print("correct")


# In[29]:


physical_states=[]
for list1 in val_master_list:
     check= all(np.abs(np.real(x)-1)<1.e-8 for x in list1)
     if check==True:
        physical_states.append(list1)
print(len(physical_states))


# ## Initial calculation with sites N=6 ##

# In[30]:


nsites=6
links= nsites-1
ndim= nsites+links
states_list= all_possible_states(ndim)
print(len(states_list))
gauss_ops = {i: gauss_law(ndim,i)for i in range(0,ndim-2,2)}
val_master_list=[]
#physical_states=[]
for state in states_list:
    val_list= np.zeros(len(gauss_ops),dtype=complex)
    for i, Gi in gauss_ops.items():
        val = (state.dag()*Gi*state)
        index=(i//2)
        val_list[index]=val
    val_master_list.append((val_list))
    state_count=[]
    for list1 in val_master_list:
     check= all(np.abs(np.real(x)-1)<1.e-8 for x in list1)
     if check==True:
        state_count.append(list1)
print(len(state_count))


# ## Physical states which satisfy Gauss Law ## 

# In[31]:


physical_states1 = []   
for state in states_list:
    val_list = np.zeros(len(gauss_ops), dtype=complex)

    for i, Gi in gauss_ops.items():
        val = (state.dag() * Gi * state)
        index = i // 2
        val_list[index] = val

    check = all(np.abs(np.real(x) - 1) < 1e-8 for x in val_list)

    if check:
        physical_states1.append(state)

print(len(physical_states1))


# In[32]:


physical_states1[-1].shape


# In[33]:


nsites=6
links=nsites-1
ndim= nsites+links
op1= total_numop(ndim)
print(op1)


# ## Physical states which are gauge invariant and have half-filled ground state ##

# In[34]:


total_numop_list=[]
for states in physical_states1 :
  numop_exp= states.dag()*total_numop(ndim)*states
  print(numop_exp)
  if(numop_exp==(nsites/2)):
    total_numop_list.append(states)
print(len(total_numop_list))


# In[35]:


total_numop_list=[]
for states in physical_states1 :
  numop_exp= states.dag()*total_numop(ndim)*states
  #print(numop_exp)
  if(numop_exp==(nsites/2)):
    total_numop_list.append(states)
print(len(total_numop_list))


# In[36]:


total_numop_list[-1]


# ## Calculation with only gauge-invariant subspace ##

# In[37]:


H_check1=  Z2_Ham_nh_new1(6,0.36,1,2*np.sqrt(0.36),0)
psi0_check = strong_coupling_vacuum(6)
el_11= physical_states1[-1].dag()*H_check1*physical_states1[-1]
print(el_11)


# In[38]:


H_phy= np.zeros((len(physical_states1),len(physical_states1)),dtype=complex)
psi_phy= np.zeros(len(physical_states1),dtype=complex)
print(H_phy.shape)
for i in range(H_phy.shape[0]):
 for j in range(H_phy.shape[0]):
    H_phy[i,j]= physical_states1[i].dag()*H_check1*physical_states1[j]

for l in range(H_phy.shape[0]):
    psi_phy[l]= physical_states1[l].dag()*psi0_check

print(H_phy)
print(psi_phy)


# In[39]:


H_phy1=    Qobj(H_phy)
psi_phy1 = Qobj(psi_phy)
tvals=  np.arange(0,10.1,0.1)
#f=1.0
#eops= [Electric_term(f,dim)]
result_check1= sesolve(H_phy1,psi_phy1,tvals)


# In[40]:


H_phy1.eigenenergies()


# In[41]:


states_check= result_check1.states
print(len(states_check))


# In[42]:


Pt_gauge_invariant=[]
for state in states_check:
    overlap= psi_phy1.dag()*state
    prob= np.linalg.norm(overlap)
    print(prob**2)
    Pt_gauge_invariant.append(prob**2)


# ## Calculation with gauge-law invariant states that are half-filled ##

# In[43]:


H_half= np.zeros((len(total_numop_list),len(total_numop_list)),dtype=complex)
psi_half= np.zeros(len(total_numop_list),dtype=complex)
print(H_half.shape)
for i in range(H_half.shape[0]):
 for j in range(H_half.shape[0]):
    H_half[i,j]= total_numop_list[i].dag()*H_check1*total_numop_list[j]

for l in range(H_half.shape[0]):
    psi_half[l]= total_numop_list[l].dag()*psi0_check

print(H_half)
print(psi_half)


# In[44]:


H_half1=    Qobj(H_half)
psi_half1 = Qobj(psi_half)
tvals=  np.arange(0,10.1,0.1)
#f=1.0
#eops= [Electric_term(f,dim)]
result_check_half1 = sesolve(H_half1,psi_half1,tvals)


# In[45]:


print(H_half1)


# In[46]:


H_half1.eigenenergies()


# In[47]:


len(H_half1.eigenstates()[1])


# In[48]:


H_half1.eigenstates()


# In[54]:


gs= (H_half1.eigenstates())[1][0]
gs_en= (H_half1.eigenenergies())[0]
eig_check1= gs_en * gs 
eig_check2= H_half1*gs
print(eig_check1 - eig_check2)


# ## Diagonal operators acting on the ground state ##

# In[56]:


def Z2_Ham_nh_check1(nsites,X,f,mu,gamma):
 vertex= nsites 
 links= nsites-1
 ndim= vertex+links
 b=(ndim/2)
 ## Constructing the Hamiltonian
 H_J=0
 H_f=0
 H_m=0
 for i in range(0,ndim-2,2):
  t1= pauliop(sigmap(),ndim,i)*0.5*pauliop(sigmax(),ndim,i+1)*pauliop(sigmam(),ndim,i+2)
  H_J+= (X)*(t1+t1.dag())
  #if 0<= i < b-2 :
   #H_J+= (X-(1j)*gamma)*(t1+t1.dag())
  #else:
    #H_J+= (X)*(t1+t1.dag())
 for j in range(1,ndim-1,2):
  H_f+= (f-(1j*gamma))*0.5*pauliop(sigmaz(),ndim,j)
 for k in range(0,ndim,2):
  vertex_index = (k) // 2
  H_m += mu * ((-1)**vertex_index)*0.5*pauliop(sigmaz(),ndim,k)

 Full_Ham= H_J+ H_f + H_m
 return(Full_Ham,H_J,H_f,H_m)


# In[64]:


H_I = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[1]
H_E = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[2]
H_m = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[3]
print(H_I)
print(H_E)
print(H_m)


# In[73]:


HI_half= np.zeros((len(total_numop_list),len(total_numop_list)),dtype=complex)
HE_half= np.zeros((len(total_numop_list),len(total_numop_list)),dtype=complex)
Hm_half= np.zeros((len(total_numop_list),len(total_numop_list)),dtype=complex)
#psi_half= np.zeros(len(total_numop_list),dtype=complex)
print(HE_half.shape)
for i in range(HE_half.shape[0]):
 for j in range(HE_half.shape[0]):
    HI_half[i,j]= total_numop_list[i].dag()*H_I*total_numop_list[j]
    HE_half[i,j]= total_numop_list[i].dag()*H_E*total_numop_list[j]
    Hm_half[i,j]= total_numop_list[i].dag()*H_m*total_numop_list[j]

print(HE_half)
print(Hm_half)


# In[74]:


#states.dag()*total_numop(ndim)*states
HI_exp1 = gs.dag()*Qobj(HI_half)*gs
HE_exp2 = gs.dag()*Qobj(HE_half)*gs
Hm_exp3 = gs.dag()*Qobj(Hm_half)*gs
print(HI_exp1,HE_exp2,Hm_exp3)
print(HI_exp1+HE_exp2+Hm_exp3)


# In[49]:


states_half_check= result_check_half1.states
print(len(states_half_check))


# In[50]:


Pt_half_filled=[]
for state in states_half_check:
    overlap1= psi_half1.dag()*state
    prob= np.linalg.norm(overlap1)
    print(prob**2)
    Pt_half_filled.append(prob**2)


# ## Plotting comparison between states which are only gauge_invariant and half-filled ground state ##

# In[112]:


plt.rcParams['xtick.labelsize'] = 22  # length of x-axis major ticks
plt.rcParams['ytick.labelsize'] = 22  # length of y-axis major ticks


# In[114]:


plt.figure(figsize=(12,8))
plt.plot(tvals,Pt_full_evolution,marker=".",label=r"$Full$ $evolution$(2048)")
plt.plot(tvals,Pt_gauge_invariant,marker=".",label=r"$Gauge$ $invariant$(64)")
plt.plot(tvals,Pt_half_filled,marker=".",label=r"$Only$ $half$ $filled$(20)")

plt.xlabel("$t$",fontsize=22)
plt.ylabel(r"$P(t)$",fontsize=22)
plt.title("$\gamma=0,Hermitian$",fontsize=22)
plt.legend()
plt.show()


# ## Diagonal operators acting on the ground state (Full Hilbert space) ##

# In[77]:


H_tot = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[0]
H_I   = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[1]
H_E   = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[2]
H_m   = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[3]

gs1= (2*H_tot.eigenstates())[1][0]
gs1_energy= 2*H_tot.eigenenergies()[0]
## expectation
HI_tot_exp = gs1.dag()*2*H_I*gs1
HE_tot_exp = gs1.dag()*2*H_E*gs1
Hm_tot_exp = gs1.dag()*2*H_m*gs1
print(HI_tot_exp,HE_tot_exp,Hm_tot_exp)
print(HI_tot_exp+HE_tot_exp+Hm_tot_exp)
print(gs1_energy)


# In[78]:


H_tot = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[0]
H_I   = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[1]
H_E   = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[2]
H_m   = Z2_Ham_nh_check1(6,0.36,1,1.2,0)[3]

gs1= (H_tot.eigenstates())[1][0]
gs1_energy= H_tot.eigenenergies()[0]
## expectation
HI_tot_exp = gs1.dag()*H_I*gs1
HE_tot_exp = gs1.dag()*H_E*gs1
Hm_tot_exp = gs1.dag()*H_m*gs1
print(HI_tot_exp,HE_tot_exp,Hm_tot_exp)
print(HI_tot_exp+HE_tot_exp+Hm_tot_exp)
print(gs1_energy)


# ## Time evolution using eigen-values and eigen-vectors ## 

# In[127]:


eig_energy_hf= H_half1.eigenstates()[0]
eig_state_hf = H_half1.eigenstates()[1]

cn_half= np.zeros(len(eig_state_hf),dtype=complex)
for l in range(len(eig_state_hf)) :
    cn_half[l] = eig_state_hf[l].dag()* psi_half1
print(cn_half)


# In[128]:


t=0
psi_state=0
for l in range(len(eig_state_hf)):
  term= cn_half[l]*np.exp(-1j*eig_energy_hf[l]*t)*eig_state_hf[l]
  psi_state = psi_state + term
print(psi_state)


# In[130]:


tvals=  np.arange(0,10.1,0.1)
Pt_eigenbasis=[]
for t in tvals:
  psi_state_t=0
  for l in range(len(eig_state_hf)):
   term= cn_half[l]*np.exp(-1j*eig_energy_hf[l]*t)*eig_state_hf[l]
   psi_state_t = psi_state_t + term
  overlap1= psi_half1.dag()*psi_state_t
  prob= np.linalg.norm(overlap1)
  print(prob**2)
  Pt_eigenbasis.append(prob**2)


# In[132]:


plt.figure(figsize=(12,8))
#Pt_eigenbasis
plt.plot(tvals,Pt_half_filled,marker=".",label=r"$Qutip$ $sesolve$")
plt.plot(tvals,Pt_eigenbasis,marker=".",label=r"$Eigen$ $value$ $method$")
plt.xlabel("$t$",fontsize=22)
plt.ylabel(r"$P(t)$",fontsize=22)
plt.title("$\gamma=0,Hermitian$",fontsize=22)
plt.legend()
plt.show()


# ## Non Hermitian expectation ##

# In[22]:


nsites= 6
dim=(2*nsites-1)
H_nh1= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0.4)
psi0 = strong_coupling_vacuum(nsites)
tvals2=  np.arange(0,10.1,0.1)
f=1.0
eops= [Electric_term(f,dim)]
result_nh= sesolve(H_nh1,psi0,tvals2,eops)
print(result_nh.expect[0])


# In[23]:


check_nh1= np.loadtxt("data_nh_1storder.dat")
check_nh2= np.loadtxt("data_nh_2ndorder.dat")
check_nh3= np.loadtxt("data_nh_TDVP.dat")
## x and y
x_checknh1= check_nh1[:,0]
y_checknh1= check_nh1[:,1]

x_checknh2= check_nh2[:,0]
y_checknh2= check_nh2[:,1]

x_checknh3= check_nh3[:,0]
y_checknh3= check_nh3[:,1]


# In[24]:


print(len(x_checknh1),len(y_checknh1))
print(len(x_checknh2),len(y_checknh2))
print(len(x_checknh3),len(y_checknh3))
print(len(tvals2),len(result_nh.expect[0]))


# In[25]:


## Plotting 
plt.figure(figsize=(12,6))
plt.plot(x_checknh1,y_checknh1,label=r"$TEBD$ $1st$ $order$")
plt.plot(x_checknh2,y_checknh2,label=r"$TEBD$ $2nd$ $order$")
plt.plot(x_checknh3,y_checknh3,label=r"$TDVP$")
plt.plot(tvals2,result_nh.expect[0],label="$Exact$")
plt.xlabel("$t$",fontsize=15)
plt.ylabel(r"$\langle E \rangle (t)$",fontsize=15)
plt.title("$\gamma=0.4,non-Hermitian")
plt.legend()
plt.show()


# ## Checking with right and left eigenvector ##

# In[86]:


nsites= 6
dim=(2*nsites-1)
H_nh1= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0.4)
psi0 = strong_coupling_vacuum(nsites)
tvals2=  np.arange(0,10.1,0.1)
f=1.0
eops= [Electric_term(f,dim)]
result_nh_r= sesolve(H_nh1,psi0,tvals2)
result_nh_l= sesolve(H_nh1.dag(),psi0,tvals2)
## statres
states_r = result_nh_r.states
states_l = result_nh_l.states
for i in range(len(states_r)):
 print(states_r[i].norm(),states_l[i].norm(),states_l[i].dag()*states_r[i])
#print(states_l[i].norm())


# In[87]:


from qutip import propagator
H_nh2= Z2_Ham_nh_new1(4,0.36,1,2*np.sqrt(0.36),0.4)
psi0_2 = strong_coupling_vacuum(4)
U  = propagator(H_nh2, tvals2)
#Ud = propagator(H_nh2.dag(), tvals2)

psiR = [U_t * psi0_2 for U_t in U]
psiL = [psi0_2.dag() * U_t.dag() for U_t in U]
for i in range(len(psiR)):
   psiR_new= psiR[i]/psiR[i].norm()
   psiL_new= psiL[i]/psiL[i].norm()
   print(psiR_new.norm(),psiL_new.norm(),psiL_new*psiR_new)
   #print(psiL_new.norm())


# In[45]:


nsites= 4
dim=(2*nsites-1)
H_nh1= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0.4)
psi0 = strong_coupling_vacuum(nsites)
tvals=  np.arange(0,10.1,0.1)
f=1.0

U  = propagator(H_nh1, tvals)
#Ud = propagator(H_nh1.dag(),tvals)

psiR = [U_t * psi0 for U_t in U]
psiL = [psi0.dag() * U_t.dag() for U_t in U]
op_exp_list=[]
op_exp_list1=[]
for i in range(len(psiR)):
   psiR_new= psiR[i]/psiR[i].norm()
   psiL_new= psiL[i]/psiL[i].norm()
   op_exp= psiL_new*Electric_term(f,dim)*psiR_new
   op_exp1= psiL[i]*Electric_term(f,dim)*psiR[i]
   print(op_exp,op_exp1)
   op_exp_list.append(op_exp)
   op_exp_list1.append(op_exp1)


# In[43]:


nsites= 4
dim=(2*nsites-1)
H_nh1= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0.4)
psi0 = strong_coupling_vacuum(nsites)
tvals2=  np.arange(0,10.1,0.1)
f=1.0
eops= [Electric_term(f,dim)]
opts = Options(normalize_output=False)

resultnh_check1= sesolve(H_nh1,psi0,tvals2,eops,options=opts)
resultnh_check2= sesolve(H_nh1,psi0,tvals2,eops)
print(resultnh_check1.expect[0])
print(resultnh_check2.expect[0])


# ## Plotting ##

# In[37]:


plt.figure(figsize=(12,6))
plt.plot(tvals2,np.real(resultnh_check2.expect[0]),label="right_exp_value")
plt.plot(tvals,op_exp_list,label="left_right_exp_value")
plt.xlabel("$t$",fontsize=15)
plt.ylabel(r"$\langle E \rangle$",fontsize=15)
plt.legend()
plt.show()


# In[46]:


plt.figure(figsize=(12,6))
plt.plot(tvals2,np.real(resultnh_check1.expect[0]),label="right_exp_value")
plt.plot(tvals,op_exp_list1,label="left_right_exp_value")
plt.xlabel("$t$",fontsize=15)
plt.ylabel(r"$\langle E \rangle$",fontsize=15)
plt.legend()
plt.show()


# In[ ]:


[-1.5        -1.49907517 -1.49654822 -1.49289505 -1.48867695 -1.48447118
 -1.48080553 -1.47810352 -1.47664495 -1.4765447  -1.47775133 -1.4800646
 -1.48316886 -1.48667791 -1.49018508 -1.49331241 -1.49575265 -1.49729977
 -1.4978649  -1.49747731 -1.49627162 -1.494464   -1.49232117 -1.49012632
 -1.48814592 -1.48660098 -1.48564552 -1.48535419 -1.48571947 -1.48665865
 -1.48802877 -1.48964746 -1.49131693 -1.49284798 -1.49408124 -1.49490332
 -1.49525633 -1.4951402  -1.49460826 -1.49375687 -1.49271119 -1.49160865
 -1.49058244 -1.48974672 -1.48918505 -1.48894326 -1.48902698 -1.48940407
 -1.49001105 -1.49076275 -1.49156371 -1.49231995 -1.49294964 -1.49339157
 -1.49361065 -1.49359988 -1.49337904 -1.49299032 -1.49249178 -1.49194956
 -1.49142965 -1.49099043 -1.49067661 -1.4905151  -1.49051335 -1.49065993
 -1.49092727 -1.49127604 -1.49166052 -1.49203436 -1.4923559  -1.49259264
 -1.49272423 -1.49274389 -1.49265813 -1.49248502 -1.49225121 -1.49198823
 -1.49172856 -1.49150178 -1.49133147 -1.49123295 -1.49121216 -1.49126569
 -1.49138188 -1.49154278 -1.49172667 -1.49191085 -1.49207431 -1.4922
 -1.49227647 -1.4922988  -1.49226861 -1.49219341 -1.49208526 -1.4919591
 -1.49183077 -1.49171513 -1.49162449 -1.49156735 -1.49154771]


# ##  EE calculation ##

# ## L=6 ##

# In[110]:


def Ent_entropy(mat):
 eiglist= mat.eigenenergies()
 EE=0
 for i in range(len(eiglist)):
    if(eiglist[i]>1.e-8):
     term= - (eiglist[i]*np.log(eiglist[i]))
     EE= EE+term
 return EE


# In[111]:


nsites= 6
dim=(2*nsites-1)
H_nh1= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0.4)
psi0 = strong_coupling_vacuum(nsites)
tvals2=  np.arange(0,30.1,0.1)
f=1.0
#eops= [Electric_term(f,dim)]
result_nh= sesolve(H_nh1,psi0,tvals2)
states_nh= result_nh.states
EE_list_check=[]
for states in states_nh:
    rhot= states*states.dag()
    rhot_norm= rhot/(states.dag()*states)
    rhoAt= rhot_norm.ptrace([0,1,2,3,4])
    #print(rhoAt.eigenenergies())
    EE_t= Ent_entropy(rhoAt)
    #print(rhot.tr(),rhoAt.tr())
    print(EE_t)
    EE_list_check.append(EE_t)


# In[113]:


EE_text1= np.loadtxt("EEnh_check_TEBD2ndorder.dat")
## x_ and y_1
x_EE1 = EE_text1[:,0]
y_EE1 = EE_text1[:,1]
## plotting
plt.figure(figsize=(12,6))
plt.plot(tvals2,EE_list_check,label="exact")
plt.plot(x_EE1,y_EE1,label="TEBD")
plt.xlabel("$t$",fontsize=15)
plt.ylabel(r"$S_({L/2},t)$",fontsize=15)
plt.legend()
plt.show()


# ## L=8 ##

# In[115]:


nsites= 8
ndim=(2*nsites-1)
H_nh_check= Z2_Ham_nh_new1(nsites,0.36,1,2*np.sqrt(0.36),0.4)
psi0 = strong_coupling_vacuum(nsites)
tvals_check=  np.arange(0,30.1,0.1)
f=1.0


# ## Gauge invariant basis ##

# In[116]:


states_list_check = all_possible_states(ndim)
print(len(states_list_check))
gauss_ops_check = {i: gauss_law(ndim,i)for i in range(0,ndim-2,2)}
val_master_list=[]


# In[ ]:




