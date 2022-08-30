from random import randrange
import numpy as np
import math

def negate(a):                      #multiplies complex by -1
    return (-a[0],-a[1])

def invert(a):                      #multiplicative inverse
    return (a[0]/(a[0]**2+a[1]**2),-a[1]/(a[0]**2+a[1]**2))

def conjugate(a):                   #returns conjugate of a complex number
    return (a[0],-a[1])

def multiplication(a,b):            #multiplies two complex number
    return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])

def division(a,b):                  #returns a/b complex
    return multiplication(a,invert(b))

def sum_cn(a,b):                    #returns sum of two numbers
    return (a[0]+a[1],b[0]+b[1])

def diff_cn(a,b):                   #returns sum of two numbers
    return (a[0]-a[1],b[0]-b[1])


def rect_to_polar(a):
    angle = math.atan(a[1]/a[0])    
    magnitude = a[0]**2 + a[1]**2
    return (magnitude,angle)

def Gauss_Schidel(V,i):  # i=1 for bus 2
    init = (0,0)
    for j in range(3):
        if j!=i:
            init = sum_cn(init,multiplication(ad_mat[i][j],Voltages_iterations[j][-1]))
            
    load_value =(0,0)
    if i==1:
        load_value = Power_2[-1]
    else:
        load_value = Power_3[-1]

    ans = multiplication(invert(ad_mat[i][i]),diff_cn(division(conjugate(load_value),conjugate(V)),init)) 

    return ans




def Q_calculator(i):
    init=(0,0)
    for j in range(3):
        init = sum_cn(init,multiplication(ad_mat[i][j],Voltages_iterations[j][-1]))

    temp = multiplication(init,conjugate(Voltages_iterations[i][-1]))
    
    Power_3.append((Power_3[-1][0],-temp[1]))

def line_loss(i,j):
    return multiplication(conjugate(ad_mat[i][j]),multiplication(Voltages_iterations[i][-1],diff_cn(conjugate(Voltages_iterations[i][-1]),conjugate(Voltages_iterations[j][-1]))))

# Parsing Transmission Line Data
TL_Data = open("rx.csv")
rx = np.genfromtxt(TL_Data, delimiter=",")
rx_new = [[(0,0) for i in range(3)] for i in range(3)]
TL_Data.close()

# Parsing Generator Data
Gen_Data = open("gen_data.csv")
Generator = np.genfromtxt(Gen_Data, delimiter=",")
Gen_new = [-1 for i in range(3)]
Gen_Data.close()

#Parsing the load Data
Load_Data = open("load_data.csv")
Load = np.genfromtxt(Load_Data, delimiter=",")
Load_new = [(0,0) for i in range(3)]
Load_Data.close()



#getting the admittance matrix
for i in range(len(rx)):
    rx_new[int(rx[i][0])-1][int(rx[i][1])-1] = (rx[i][3],rx[i][4])
    rx_new[int(rx[i][1])-1][int(rx[i][0])-1] = (rx[i][3],rx[i][4])

ad_mat = [[(0,0) for i in range(3)] for i in range(3)]

for i in range(3):
    for j in range(3):
        if i!=j:
            ad_mat[i][j] = negate(invert(rx_new[i][j]))

for i in range(3):
    a = -sum([ad_mat[i][j][0] for j in range(3) if i!=j])
    b = -sum([ad_mat[i][j][1] for j in range(3) if i!=j])
    ad_mat[i][i] = (a,b)

for i in ad_mat:
    print(*i)

# Getting the Generator Data
for i in range(len(Generator)):
    Gen_new[int(Generator[i][0])-1] = Generator[i][2]

# getting the Load data
for i in range(3):
    Load_new[i] = (Load[i][2],Load[i][3])

# initial guess of Voltages
Voltages_iterations = [[(1.05,0)],
                       [(1,0)],
                       [(1.04,0)]]


Power_2 = [(4,2)]
Power_3 = [(2,-1)]

iterations = 0

condition = True
tolerance = 1e-4

while condition:
    for j in range(3):
        if j==0:
            # Slack bus The voltage value remains same
            temp = Voltages_iterations[j][-1] 
            Voltages_iterations[j].append(temp)           
        elif j==1:
            # PQ bus
            temp = Power_2[-1]
            Power_2.append(temp)

            return_value = Gauss_Schidel(Voltages_iterations[j][-1],j)
            Voltages_iterations[j].append(return_value)

        else:
            # PV bus
            Q_calculator(j)
            return_value = Gauss_Schidel(Voltages_iterations[j][-1],j)

            angle = math.atan(return_value[1]/return_value[0])

            new_value = (Voltages_iterations[j][0][0]*math.cos(angle),Voltages_iterations[j][0][0]*math.sin(angle))

            Voltages_iterations[j].append(new_value)
    
    iterations -= -1

    if abs(rect_to_polar(Voltages_iterations[1][-1])[0] - rect_to_polar(Voltages_iterations[1][-2])[0])<=tolerance:
        condition=False

print()

for i in range(3):
    for j in range(3):
        if i>j:
            print("The line loss of line",i,j,"is", line_loss(i,j))


print()

print("The model converged after",iterations,"iterations.")

print()

print("For the PQ bus:")
print("The magnitude of Voltage is",rect_to_polar(Voltages_iterations[1][-1])[0])
print("The phase of Voltage is",rect_to_polar(Voltages_iterations[1][-1])[1],"radians.")

print()

print("For the PV bus:")
print("The magnitude of Voltage is",1.04,".")
print("The phase of Voltage is",rect_to_polar(Voltages_iterations[2][-1])[1],"radians.")
