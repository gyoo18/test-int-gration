import sympy
def detA(i,j,N):
    X=[-2,0,5,12]
    Y=[-22,2,27]
    if i==N-1:
        return X[i]^j
   
   
    somme=0
    for k in range(j):
        print(i,j,N,k)
        somme+=X[i+1]**(k-1)*(-1)**(1+k)*detA(i+1,k,N)
       
    for k in range(j+1,N):
        print(i,j,N,k)
        somme+=X[i+1]**(k-1)*(-1)**(1+k)*detA(i+1,k,N)
    return somme
x1=sympy.symbols('x1')
x2=sympy.symbols('x2')
x3=sympy.symbols('x3')
x4=sympy.symbols('x4')
D=(x2*x3**2-x2**2*x3)-x1*(x3**2-x2**2)+x1**2*(x3-x2)
Df=-(x1-x2)*(x1-x3)*(x2-x3)
P4=(x1-x2)*(x1-x3)*(x1-x4)*(x2-x3)*(x2-x4)*(x3-x4)
P2=(1*x2-x1*1)
P2=-(x1-x2)
P1=1
def factoriel(N):
    if N==0:
        return 1
    return N*factoriel(N-1)

def DetX(X,N):
    if len(X)==1:
        return 1

    produ=1
    for i in range(1,len(X)):
        produ*=(X[0]-X[i])
       
    return produ*DetX(X[1:],N)

def MaP(X):
    Ma=[]
    N=len(X)
    for i in range(N):
        row=[]
        for j in range(N):
            row.append(X[i]**(j))
        Ma.append(row)
    Ma=sympy.Matrix(Ma)
    return Ma

def MYP(X,Y,c):
    Ma=[]
    N=len(X)
    for i in range(N):
        row=[]
        for j in range(N):
            if j==c:
                row.append(Y[i])
            else:
                row.append(X[i]**(j))
        Ma.append(row)
    Ma=sympy.Matrix(Ma)
    return Ma
ma=MaP([-2,0,5,12])

def Poly(X,Y):
    a=[]
    d=MaP(X).det()
    for c in range(len(Y)):
        a.append(MYP(X,Y,c).det()/d)
    return a
