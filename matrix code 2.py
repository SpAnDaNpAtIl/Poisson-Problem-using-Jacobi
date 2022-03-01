import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

start = time.time()
h = 1/10 #change the h for different distribution of points
N = int(1/h)
f=0
A = np.zeros([(N-1)**2, (N-1)**2])
u = np.zeros((N-1)**2)
B= np.zeros((N-1)**2)
ucombiner=[]
ucombiner.append(u)
T = np.zeros([N-1, N-1])
I = np.zeros([N-1, N-1])
for i in range(N-1):
    for j in range(N-1):
        if(i==j):
            T[i][j] = 4
            try:
                T[i][j-1] = -1
            except:
                None
            try:
                T[i][j+1] = -1
            except:
                None
                
    T[0][N-2]=0 #because in python arr[-1] is the last element and not element outside the array
                
for i in range(N-1):
    for j in range(N-1):
        if(i==j):
            I[i][i] = -1
            
for i in range(N-1):
    A[(N-1)*i:(N-1)*i+T.shape[0], (N-1)*i:(N-1)*i+T.shape[1]] += T
    try:
        A[(N-1)*i:(N-1)*i+T.shape[0], (N-1)*(i+1):(N-1)*(i+1)+T.shape[1]] += I
    except:
        None
    if(i>0):
        try:
            A[(N-1)*i:(N-1)*i+T.shape[0], (N-1)*(i-1):(N-1)*(i-1)+T.shape[1]] += I
        except:
            None
        
#matrix done similar using T and -I
#defining u00 vector containing index data of u vector(like (1,2) (3,4))
u00 = []
for i in range(N-1):
    for j in range(N-1):
        u00.append([j+1, i+1])
        
        
for i in range(len(u00)):
    if(u00[i][0]==1):
        None #since g is 0 for x=0     
    if(u00[i][0]==N-1):
        B[i] += h*u00[i][1] #x=Nth i.e. right side closer point      
    if(u00[i][1]==1):
        B[i] += (h*u00[i][0] - 1)*np.sin(h*u00[i][0])  #y=0 closer points  
    if(u00[i][1]==N-1):
        B[i] += (h*u00[i][0])*(2 - h*u00[i][0]) #y=nth upper closer points

        
#print(u00, B)
#Solving the u's for each iteration
def productor(lis1, lis2):
    return sum(i[0] * i[1] for i in zip(lis1, lis2))

def ratioer(lis1, lis2):
    num = sum((i[0] - i[1])**2 for i in zip(lis1, lis2))
    den = sum(i**2 for i in lis1)
    return (num/den)**0.5

x = []
y=[]
a=1
booler = True #booler will not switch the value till relative error is less than 10^-8
while(booler == True):
    unew = np.zeros((N-1)**2)
    x.append(a)
    for i in range(len(unew)):
        unew[i] += (B[i] - productor(A[i], ucombiner[-1]) + ucombiner[-1][i]*A[i][i])/A[i][i]
    y.append(math.log10(ratioer(unew, ucombiner[-1]))) 
    if(ratioer(unew, ucombiner[-1])< 10**(-8)):
        booler = False   
    ucombiner.append(unew)
    print(a)
    if(len(ucombiner)>2):
        ucombiner.pop(0) #removing unwanted data to save memory
    a+=1
    
    
#print(ucombiner[-1])
plt.plot(x,y)
plt.xlabel('No of iterations')
plt.ylabel('relative error')
plt.legend(['h='+str(h)])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x1 = [h*(i+1) for i in range(N-1)]*(N-1)
y1=[]
for i in range(N-1):
    y1 += [i+1]*(N-1)
    
# print(y1)
#print(ucombiner[-1])
ax.plot_trisurf(x1, y1, ucombiner[-1],
                linewidth = 0.2,
                antialiased = True);

plt.xlabel('x')
plt.ylabel('y')
ax.set_zlabel('u')
plt.show()
end = time.time()
print('Time required for computation is '+str(end-start)+'seconds and no. of total iterations were '+str(a-1))
#h=1/40 takes a lot of time to compute. If iteration number is printed in the code, we can see that part of code till matrix generation is pretty fast.
#But finding individual u_i,j takes a lot of time. 

    
    
    

        
    
