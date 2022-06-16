import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
import csv
from alive_progress import alive_bar
#################################################################################
#defining what air is
air = "O2:0.21,N2:0.79"
#simulation options
Duration = 250
iter1=200 #Quality of products of combustion calculations
iter2=500 #Quality of pressure and temperature calculations
plot_gen = True
Outside_data = ''
#Idle variables (placeholders)
low_limit = 0
high_limit = 1
low_limit_fig=0
high_limit_fig=1
ox_R = 280
fuel_R = 280

#Mechanism used for the process
gas = ct.Solution('gri30.yaml')
#Gas state
P_0 = ct.one_atm
T_0 = 1000
gas.TP = T_0, P_0
#Choosing oxidizer and fuel
fuel = input('Fuel:\n air and O2 - H2, CO, CH4, C2H6, C3H8\n only air - CH3OH\n Choose: ')
if fuel == 'CH3OH':
    oxidizer = 'air'
else:
    oxidizer = input('Oxidizer (air or O2): ')

#Loading flammability limits
if oxidizer == 'O2':
    Outside_data = 'O2.csv'
    ox_R = 260
elif oxidizer == 'air':
    Outside_data = 'AIR.csv'
    oxidizer = air
    ox_R = 287
with open(Outside_data) as f:
    reader = csv.reader(f, delimiter='\t')
    Limits = [(col1, float(col2), float(col3),float(col4))
                for col1, col2, col3, col4 in reader]
for a in range(len(Limits)):
        arr_limits=np.asarray(Limits)
        if arr_limits[a,0] == fuel :
            low_limit= float(arr_limits[a,1])
            high_limit = float(arr_limits[a, 2])
            fuel_R = float(arr_limits[a, 3])

#converting volume fractions to mass fractions
low_limit_fig= (low_limit/fuel_R) /(((100-low_limit)/ox_R)+(low_limit/fuel_R))
high_limit_fig= (high_limit/fuel_R) /(((100-high_limit)/ox_R)+(high_limit/fuel_R))
low_limit = 0
high_limit =1

#Allocating space for data
fo = 0 #Fuel to all Mixture ratio
fos=np.zeros(iter1)
data = np.zeros((iter1,10))
data1 = np.zeros((iter1,2))
arr3 = np.zeros(10)
#Main working loop
with alive_bar(iter1,force_tty=True) as bar:
    for m in range(iter1):
        fo += (high_limit-low_limit)/iter1
        if fo >=high_limit:
            fo=high_limit-0.0001
        #Main calculations
        gas.TP = T_0, P_0
        gas.set_mixture_fraction(fo, fuel, oxidizer)
        r = ct.IdealGasReactor(gas)
        sim = ct.ReactorNet([r])
        time = Duration/iter2
        times = np.zeros(iter2)

        #Getting max T and max P
        p_max = P_0
        T_max = T_0
        for n in range(iter2):
            time += Duration/iter2
            sim.advance(time)
            times[n] = time # time in s
            if T_max < r.T:
                T_max = r.T
            if p_max < r.thermo.P:
                p_max = r.thermo.P
# Bubble Sort for 10 most present substances
        if m ==0:
            gas.set_equivalence_ratio(1.0, fuel, oxidizer)
            name = gas.species_names
            string2 = gas.Y
            arr1 = np.array(name)
            arr2 = np.array(string2)
            n = len(arr2)
            swapped = False
            for i in range(n - 1):
                for j in range(0, n - i - 1):
                    if arr2[j] < arr2[j + 1]:
                        swapped = True
                        arr2[j], arr2[j + 1] = arr2[j + 1], arr2[j]
                        arr1[j], arr1[j + 1] = arr1[j + 1], arr1[j]
            fo=low_limit
            fos[m] = fo
            arr3=arr1[:10]
        else:
             fos[m] = fo
             data[m,0:] = r.thermo[arr3].Y
             data1[m,0] = T_max
             data1[m,1] = p_max
        bar()

# Plot the results if matplotlib is installed.
if plot_gen==True:
    plt.xlabel('Mass of fuel in the mixture before combustion[%]')
    plt.ylabel('Products of combustion [%]')
    for o in range(10):
        plt.plot(fos, data[:,o], label=arr3[o])
    plt.legend()
    plt.show()

    plt.clf()
    plt.subplot(2,1, 1)
    plt.plot(fos, data1[:, 0], label='Calculated conditions')
    plt.plot([0,1], [T_0,T_0], label='Starting conditions')
    #plt.plot([low_limit_fig, low_limit_fig], [T_0,T_max], color='C3')
    #plt.plot([high_limit_fig, high_limit_fig], [T_0,T_max], color='C3')
    plt.xlabel('Mass of fuel in the mixture [%]')
    plt.ylabel('Maximum temperature (K)')
    plt.legend()
    plt.subplot(2,1, 2)
    plt.plot(fos, data1[:, 1]/100000, label='Calculated conditions')
    plt.plot([0,1], [P_0/100000,P_0/100000],label='Starting conditions')
    #plt.plot([low_limit_fig, low_limit_fig], [P_0/100000,p_max/100000], color='C3')
    #plt.plot([high_limit_fig, high_limit_fig], [P_0/100000,p_max/100000], color='C3')
    plt.xlabel('Mass of fuel in the mixture [%]')
    plt.ylabel('Maximum pressure (bar)')
    plt.legend()
    plt.show()







