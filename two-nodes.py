import numpy as np
import math as m
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import minimize
plt.rcdefaults()

Boltz=5.67*10**(-8)     #Cte de Boltzmann
S_ext=0.06
a=0.34                  #Valor promedio del albedo
T_p=288                 #Temperatura planetaria promedio

E_s=0.9                 #Emisividad de la superficie
E_i=0.05                #Emisividad interna

F_a=0.5                 #Factor de vista para el cálculo del albedo
F_s=0.5                 #Factor de vista entre ??? y el planeta

K_FR4=0.25              #Conductividad calorífica del FR-4

A_s=0.03231             #Área que encara al sol
A_p=0.03231             #Área que encara al planeta (albedo e IR)
A_r=0.06462             #Área radiativa
alpha_s=0.65            #Coeficiente de absortividad

Coeficientes_1=np.genfromtxt('data/Coeficientes.csv', delimiter=',')
Coeficientes_2=np.genfromtxt('data/Coeficientes2.csv', delimiter=',')

AvCf1=np.average(Coeficientes_1)
AvCf2=np.average(Coeficientes_2)

C_int=349.5046
C_ext=366.1146
q_int=0.2                #Potencia disipada por la batería en forma de calor

#Cálculos realizados al interior del modelo
r_ext=(S_ext/(4*m.pi))**(0.5)   #Radio exterior ficticio
r_int=0.01                      #Radio interior ficticio

K=(4*m.pi*K_FR4*r_int*r_ext)/(r_ext-r_int)
R=(4*m.pi*r_int**2)/((1/E_i)+((1-E_s)/(E_s))*(r_int/r_ext)**2)

#Mitad del período radiativo
TAN=2820

#Resolución de la ecuación diferencial involucrada
def TheSuchaiSystem(state,t):
    T_int, T_ext = state
    q_int=0.2
    if t%(TAN*2)<TAN:
        J_s=1391.5
    elif t%(TAN*2)>TAN: 
        J_s=0
    
    if t>(138260) and t<(144260):
        q_int=0.8
    elif t>=(144260):
        q_int=0.2

    q_s=A_s*alpha_s*J_s
    q_a=A_p*F_a*alpha_s*J_s*a
    q_p=A_p*F_s*E_s*Boltz*(T_p**4)
    q_r=A_r*E_s*Boltz*(T_ext**4)
    dT_int=((K*(T_ext-T_int)+R*Boltz*((T_ext**(4))-(T_int**(4)))+q_int)/C_int)
    dT_ext=((K*(T_int-T_ext)+R*Boltz*((T_int**(4))-(T_ext**(4)))+q_s+q_a-q_r+q_p)/C_ext)
    return [dT_int,dT_ext]

t=np.arange(0,200000,60)   #El vector del tiempo

cin=[273.15,273.15]   #Condiciones iniciales del sistema de ecuaciones diferenciales

sol_1=odeint(TheSuchaiSystem,cin,t)

sol=sol_1-273.15

sol_pedazo=sol[int(((TAN/48)*40-25)):, 0]
np.savetxt("Vector.csv", sol_pedazo, delimiter=",")

#Sección donde se van creando todos los vectores
#Lector del CSV y de datos
datos_temp=np.genfromtxt('data/tm_status-18-06-26.csv', delimiter=',')
datos_temp_2=np.genfromtxt('data/tm_status-18-06-30.csv', delimiter=',')
datos_temp_3=np.genfromtxt('data/tm_status-18-06-16.csv', delimiter=',')
datos_temp_4=np.genfromtxt('data/tm_status-18-06-01.csv', delimiter=',')
datos_temp_1=np.genfromtxt('data/tm_status-18-01-07.csv', delimiter=',')
datos_temp_5=np.genfromtxt('data/tm_status-18-03-05.csv', delimiter=',')
datos_temp_exp=np.genfromtxt('data/tm_status-17-10-27.csv', delimiter=',')

#gyvg
Vect_cof1=C_int*np.ones(len(Coeficientes_1))
Vect_cof2=C_ext*np.ones(len(Coeficientes_2))

plt.figure()
plt.ylim((0, 600))
plt.plot(np.linspace(1,len(Coeficientes_1),len(Coeficientes_1)),Coeficientes_1,label='C_int',color='b')
plt.plot(np.linspace(1,len(Coeficientes_2),len(Coeficientes_2)),Coeficientes_2,label='C_ext',color='orange')
plt.plot(np.linspace(1,len(Coeficientes_2),len(Coeficientes_2)),Vect_cof1,label='Average C_int',color='b')
plt.plot(np.linspace(1,len(Coeficientes_2),len(Coeficientes_2)),Vect_cof2,label='Average C_ext',color='orange')
plt.xlabel("Number of Dataset")
plt.ylabel("Heat Capacity [J/K]")
plt.grid()
plt.legend(['$C_{int}$', '$C_{ext}$','$Average\; C_{int}$','$Average\; C_{ext}$'])
plt.show()

temp_np=datos_temp[2:677,1:]
temp_np_1=datos_temp_1[2:192,1:]
temp_np_2=datos_temp_2[2:417,1:]
temp_np_3=datos_temp_3[2:,1:]
temp_np_4=datos_temp_4[2:,1:]
temp_np_5=datos_temp_5[2:,1:]
temp_np_exp=datos_temp_exp[2:,1:]

#Creador de temperaturas promedio
temp_pr=np.zeros(len(temp_np))
for i in range(len(temp_pr)):
    temp_pr[i]=(temp_np[i,0]+temp_np[i,1]+temp_np[i,2]+temp_np[i,3])/4
temp_df=temp_pr[43:571]
vec_t=np.zeros(len(temp_df))
for j in range(len(vec_t)):
    vec_t[j]=1+60*j
vec_t_l=list(vec_t)
temp_df_l=list(temp_df)

del(vec_t_l[len(sol_pedazo):])
del(temp_df_l[len(sol_pedazo):])

temp_pr_1=np.zeros(len(temp_np_1))
for i in range(len(temp_pr_1)):
    temp_pr_1[i]=(temp_np_1[i,0]+temp_np_1[i,1]+temp_np_1[i,2]+temp_np_1[i,3])/4
temp_df_1=temp_pr_1[37:]
vec_t_1=np.zeros(len(temp_df_1))
for j in range(len(vec_t_1)):
    vec_t_1[j]=1+60*j
vec_t_l_1=list(vec_t_1)
temp_df_l_1=list(temp_df_1)


temp_pr_2=np.zeros(len(temp_np_2))
for i in range(len(temp_pr_2)):
    temp_pr_2[i]=(temp_np_2[i,0]+temp_np_2[i,1]+temp_np_2[i,2]+temp_np_2[i,3])/4
temp_df_2=temp_pr_2[17:]
vec_t_2=np.zeros(len(temp_df_2))
for j in range(len(vec_t_2)):
    vec_t_2[j]=1+60*j
vec_t_l_2=list(vec_t_2)
temp_df_l_2=list(temp_df_2)


temp_pr_3=np.zeros(len(temp_np_3))
for i in range(len(temp_pr_3)):
    temp_pr_3[i]=(temp_np_3[i,0]+temp_np_3[i,1]+temp_np_3[i,2]+temp_np_3[i,3])/4
temp_df_3=temp_pr_3[21:]
vec_t_3=np.zeros(len(temp_df_3))
for j in range(len(vec_t_3)):
    vec_t_3[j]=1+60*j
vec_t_l_3=list(vec_t_3)
temp_df_l_3=list(temp_df_3)


temp_pr_4=np.zeros(len(temp_np_4))
for i in range(len(temp_pr_4)):
    temp_pr_4[i]=(temp_np_4[i,0]+temp_np_4[i,1]+temp_np_4[i,2]+temp_np_4[i,3])/4
temp_df_4=temp_pr_4[21:]
vec_t_4=np.zeros(len(temp_df_4))
for j in range(len(vec_t_4)):
    vec_t_4[j]=1+60*j
vec_t_l_4=list(vec_t_4)
temp_df_l_4=list(temp_df_4)


temp_pr_5=np.zeros(len(temp_np_5))
for i in range(len(temp_pr_5)):
    temp_pr_5[i]=(temp_np_5[i,0]+temp_np_5[i,1]+temp_np_5[i,2]+temp_np_5[i,3])/4
temp_df_5=temp_pr_5[48:]
vec_t_5=np.zeros(len(temp_df_5))
for j in range(len(vec_t_5)):
    vec_t_5[j]=1+60*j
vec_t_l_5=list(vec_t_5)
temp_df_l_5=list(temp_df_5)

temp_pr_exp=np.zeros(len(temp_np_exp))
for i in range(len(temp_pr_exp)):
    temp_pr_exp[i]=(temp_np_exp[i,0]+temp_np_exp[i,1]+temp_np_exp[i,2]+temp_np_exp[i,3])/4
temp_df_exp=temp_pr_exp[:]
vec_t_exp=np.zeros(len(temp_df_exp))
for j in range(len(vec_t_exp)):
    vec_t_exp[j]=1+60*j
vec_t_l_exp=list(vec_t_exp)
temp_df_l_exp=list(temp_df_exp)

temp_df_l_exp=temp_df_l_exp[:(len(temp_df_l_exp)-3)]

plt.figure()
plt.scatter(np.linspace(1,len(sol_pedazo[:len(temp_df_l_exp)]),len(sol_pedazo[:len(temp_df_l_exp)])),sol_pedazo[:len(temp_df_l_exp)],label='Model output temperature')
plt.scatter(np.linspace(1,len(temp_df_l_exp),len(temp_df_l_exp)),temp_df_l_exp,label='Average measured temperature')
plt.plot(np.linspace(1,len(temp_np_exp[:(len(temp_np_exp)-3),2]),len(temp_np_exp[:(len(temp_np_exp)-3),2])),temp_np_exp[:(len(temp_np_exp)-3),2],color='green',label='Coldest measured temperature')
plt.plot(np.linspace(1,len(temp_np_exp[:(len(temp_np_exp)-3),3]),len(temp_np_exp[:(len(temp_np_exp)-3),3])),temp_np_exp[:(len(temp_np_exp)-3),3],color='red',label='Hottest measured temperature')
plt.xlabel("Time [min]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.scatter(np.linspace(1,len(sol_pedazo[:len(temp_df_l_1)]),len(sol_pedazo[:len(temp_df_l_1)])),sol_pedazo[:len(temp_df_l_1)],label='Model output temperature',marker='.')
plt.scatter(np.linspace(1,len(temp_df_l_exp),len(temp_df_l_exp)),temp_df_l_exp,label='Average measured temperature',marker='.')
plt.xlabel("Time [min]")
plt.ylabel("Temperature [°C]")
plt.grid()
plt.legend()
plt.show()


#Cálculo 
def RMS(C):
    C_int=C[0]
    C_ext=C[1]
    
    def TheSuchaiSystem(state,t):
        T_int, T_ext = state
        q_int=0.2
        if t%(TAN*2)<TAN:
            J_s=1391.5
        elif t%(TAN*2)>TAN: 
            J_s=0

        if t>141160 and t<147660:
            q_int=0.8
        elif t>=147660:
            q_int=0.2
        
        q_s=A_s*alpha_s*J_s
        q_a=A_p*F_a*alpha_s*J_s*a
        q_p=A_p*F_s*E_s*Boltz*(T_p**4)
        q_r=A_r*E_s*Boltz*(T_ext**4)
        dT_int=((K*(T_ext-T_int)+R*Boltz*((T_ext**(4))-(T_int**(4)))+q_int)/C_int)
        dT_ext=((K*(T_int-T_ext)+R*Boltz*((T_int**(4))-(T_ext**(4)))+q_s+q_a-q_r+q_p)/C_ext)
        return [dT_int,dT_ext]
    
    cin=[273.15,273.15]   #Condiciones iniciales del sistema de ecuaciones diferenciales
    t=np.arange(0,200000,60)    #El vector del tiempo
    sol_1=odeint(TheSuchaiSystem,cin,t)
    sol=sol_1-273.15
    sol_pedazo=sol[int(((TAN/48)*40+10))-35:, 0]       
    return sum(((sol_pedazo[:len(temp_df_l_exp)]-temp_df_l_exp[:len(sol_pedazo[:len(temp_df_l_exp)])])**2)/(len(sol_pedazo[:len(temp_df_l_exp)])))

MSE=sum((((sol_pedazo[:len(temp_df_l)]-temp_df_l[:len(sol_pedazo[:len(temp_df_l)])])**2)/(len(sol_pedazo[:len(temp_df_l)]))))
MSE_1=sum((((sol_pedazo[:len(temp_df_l_exp)]-temp_df_l_exp[:len(sol_pedazo[:len(temp_df_l_exp)])])**2)/(len(sol_pedazo[:len(temp_df_l_exp)]))))

C_0=np.array([291,401])
res=minimize(RMS,C_0, method='nelder-mead',options={'xtol': 1e-8, 'disp': True})
print(res)

