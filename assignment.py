class LPsolver:

    # Class LPsolver created by --
    # Name  : UJJWAL SINGH
    # S.no. : 70

    def __init__(self):                                                     # Initialising instance variables
        self.objective=None
        self.constraints=None
        self.table=None
        self.result=None
    
    def solve(self,L,M,option):
        return (self.Simplex(L,M))

    def Simplex(self,L,M):                                                  # Simplex function takes two lists as arguments
            
        # List L contains co-efficients of variables in the objective
        # e.g.- If objective is (2x+5y+7z) then L=[2,5,7]
        # List M contains co-efficients of variables in the constraints and the last element is the constant
        # List M is a nested list
        # e.g. - If the constraints are x<=200 , y<=10 , z<=200 and x+2y+3z<=700 then --  
        # M=[[1,0,0,200],[0,1,0,10],[0,0,1,200],[1,2,3,700]], i.e.,
        # 1.x+0.y+0.z<=200 , 0.x+1.y+0.z<=10 , 0.x+0.y+1.z<=200 and 1.x+2.y+3.z<=700
        
        self.objective=L
        self.constraints=M
        no_of_variables=len(self.objective)                                 # variable to store total number of variables
        no_of_equations=len(self.constraints)                               # variable to store total number of equations(constraints)
        self.table=[[0 for i in range(no_of_variables+no_of_equations+2)] for j in range(no_of_equations+1)]   # initialising the augmented matrix
        for i in range(no_of_equations):
            for j in range(no_of_variables):
                self.table[i][j]=self.constraints[i][j]                     # filling up the matrix
        for i in range(no_of_variables):
            self.table[no_of_equations][i]=(-1)*self.objective[i]           # filling up the matrix
        for i in range(no_of_equations+1):
            self.table[i][no_of_variables+i]=1                              # filling up the matrix
        for i in range(no_of_equations):
            self.table[i][no_of_variables+no_of_equations+1]=self.constraints[i][-1]              # filling up the matrix
        self.table[-1][-1]=0                                                # matrix ready
        for i in range(no_of_equations+1):
            self.table[i][no_of_variables+i]=self.table[i][-1]/self.table[i][no_of_variables+i]
        flag=0
        if(min(self.table[no_of_equations])>=0):                            # testing terminal condition
            flag=1
        while(flag==0):                                                     # loop
            pivot_col=self.table[no_of_equations].index(min(self.table[no_of_equations])) # finding column of pivot
            pivot_row,temp,flag2=-1,0,0                                         
            for i in range(no_of_equations):                                                        
                if(self.table[i][pivot_col]>0 and (pivot_row==-1 or self.table[i][-1]/self.table[i][pivot_col]<temp)):   # finding row of pivot
                    pivot_row=i
                    temp=self.table[i][-1]/self.table[i][pivot_col]
                    flag2=1
            if(flag2==0):
                return ("Unbound")                                          # terminal condition 
            pivot=self.table[pivot_row][pivot_col]
            for i in range(no_of_equations+1):                              # loop to set all other elements of pivot's column to 0 by row operations
                if(i==pivot_row):
                    continue
                temp=self.table[i][pivot_col]/pivot
                for j in range(no_of_variables+no_of_equations+2):
                    if(j==pivot_col):
                        self.table[i][j]=0.0
                    else:    
                        self.table[i][j]=self.table[i][j]-(temp*self.table[pivot_row][j])
            if(min(self.table[no_of_equations])>=0):
                flag=1                                                      # terminal condition
        self.result=""                                                      # string to store final answer
        for i in range(no_of_variables):                                    # building up the string
            count,temp,flag=0,0,1
            for j in range(no_of_equations+1):
                if(self.table[j][i]!=0):
                    count+=1
                    if(count>1):
                        flag=0
                        break
                    temp=self.table[j][-1]/self.table[j][i]
            if(flag==1):
                self.result+="x_"+str(i+1)+" = "+str(temp)+" ; "
            else:
                self.result+="x_"+str(i+1)+"=0 ; "
        self.result+="Maximum value of objective = "+str(self.table[-1][-1])
        return (self.result)                                                # return result

lps=LPsolver()                                                              # object creation
solution=lps.solve([1,3],[[1,2,10],[1,0,5],[0,1,4]],option="Simplex")       # calling Simplex method
print(solution)

# New Class #

class LinearSystemSolver:

    # Class LinearSystemSolver created by --
    # Name  : SUSHISH KUMAR
    # S.no. : 69

    def __init__(self):                                                     # initialising instance variables
        self.A=None
        self.b=None
        self.max_iterations=None
        self.accuracy=None
        self.initial_guess=None
        self.augmented=None
        self.result=None
    
    def solve(self,A,B,method,max_iterations=10000,accuracy=0.0001,initial_guess=None):
        
        # max_iterations, accuracy, and initial_guess are optional arguments
        
        if(method=="gauss"):
            return (self.Gauss(A,B))
        elif(method=="gauss-jordan"):
            return (self.GaussJordan(A,B))
        else:
            return (self.GaussSiedel(A,B,max_iterations,accuracy,initial_guess))

    def Gauss(self,A,B):                                                    # Gauss function takes two lists as arguments
                                                                            # A is the co-efficient matrix, B is matrix of constants

        # Suppose we have to solve the system --
        # 1.x+1.y+1.z=10
        # 0.x+1.y+0.z=2
        # 0.x+0.y+1.z=5
        # Then --
        # A=[[1,1,1],[0,1,0],[0,0,1]]
        # B=[10,2,5]
        
        self.A=A
        self.b=B
        n=len(self.A)                                                       # finding order of co-efficient matrix
        self.augmented=[[0 for i in range(n+1)] for j in range(n)]          # initialising the augmented matrix
        for i in range(n):
            for j in range(n):
                self.augmented[i][j]=self.A[i][j]                           # filling up the entries of matrix A
        for i in range(n):
            self.augmented[i][n]=self.b[i]                                  # filling up the entries of matrix b
        for i in range(n-1):
            if(self.augmented[i][i]==0):                                    # checking if the pivot is 0
                flag=0  
                for j in range(i+1,n):
                    if(self.augmented[j][i]!=0):                            # swapping rows if pivot is 0
                        self.augmented[i],self.augmented[j]=self.augmented[j],self.augmented[i]          
                        flag=1
                        break
                if(flag==0):                                                
                    continue                                                # continue if the column contains only 0s
            for j in range(i+1,n):
                temp=self.augmented[j][i]/self.augmented[i][i]                                    
                for k in range(n+1):                                        # performing row operations to make the entries below the pivot=0
                    if(k==i):
                        self.augmented[j][k]=0.0
                    else:
                        self.augmented[j][k]=self.augmented[j][k]-temp*self.augmented[i][k]
        self.result=[0 for i in range(n)]                                   # initialisng the matrix containing the solution matrix
        for i in range(n-1,-1,-1):
            s=0
            for j in range(n):                                              # solving using back substitution
                s+=self.augmented[i][j]*self.result[j]
            self.result[i]=(self.augmented[i][n]-s)/self.augmented[i][i]
        return(self.result)                                                 # returning solution matrix

    def GaussJordan(self,A,B):                                              # GaussJordan function takes two lists as arguments
                                                                            # A is the co-efficient matrix, B is matrix of constants

        # Suppose we have to solve the system --
        # 1.x+1.y+1.z=6
        # 2.x+1.y+(-1).z=1
        # (-1).x+2.y+2.z=9
        # Then --
        # A=[[1,1,1],[2,1,-1],[-1,2,2]]
        # B=[6,1,9]
        
        self.A=A
        self.b=B
        n=len(self.A)                                                       # finding order of co-efficient matrix
        self.augmented=[[0 for i in range(n+1)] for j in range(n)]          # initialising the augmented matrix
        for i in range(n):
            for j in range(n):
                self.augmented[i][j]=self.A[i][j]                           # filling up the entries of matrix A
        for i in range(n):
            self.augmented[i][n]=self.b[i]                                  # filling up the entries of matrix b
        for i in range(n):
            if(self.augmented[i][i]==0):                                    # checking if the pivot is 0
                flag=0
                for j in range(n):
                    if(self.augmented[j][i]!=0):                            # swapping rows if pivot is 0          
                        self.augmented[i],self.augmented[j]=self.augmented[j],self.augmented[i]
                        flag=1
                        break
                if(flag==0):
                    continue                                                # continue if the column contains only 0s
            for j in range(n):
                if(j==i):
                    continue
                temp=self.augmented[j][i]/self.augmented[i][i]
                for k in range(n+1):                                        # performing row operations to make the entries in the the pivot column=0
                    if(k==i):
                        self.augmented[j][k]=0.0
                    else:
                        self.augmented[j][k]=self.augmented[j][k]-temp*self.augmented[i][k]
        self.result=[0 for i in range(n)]                                   # initialisng the matrix containing the solution matrix
        for i in range(n):                                                  
            self.result[i]=self.augmented[i][n]/self.augmented[i][i]        # solving
        return(self.result)                                                 # returning solution matrix
    
    def GaussSiedel(self,A,B,max_iterations,accuracy,initial_guess):        # GaussSiedel function takes two lists 
                                                                            # (A is the co-efficient matrix, B is matrix of constants),         
                                                                            # maximmum allowed iterations, desired accuracy and initial guess as argument

        # Suppose we have to solve the system --
        # 1.x+1.y+1.z=10
        # 0.x+1.y+0.z=2
        # 0.x+0.y+1.z=5
        # Then --
        # A=[[1,1,1],[0,1,0],[0,0,1]]
        # B=[10,2,5]
        # max_iterations (optional argument) is the maximum number of iterations 
        # accuracy (optional argument) is the desired amount of Accuracy
        # initial_guess (optional argument) is the starting matrix, e.g.--
        #initial_guess=[1,1,1]
        
        from numpy import matrix                                            
        from numpy.linalg import inv
        self.A=A
        self.b=B
        self.max_iterations=max_iterations
        self.accuracy=accuracy
        if(initial_guess!=None):
            self.initial_guess=initial_guess
        else:
            self.initial_guess=[1 for i in range(len(self.A))]
        n=len(self.A)                                                       # finding order of co-efficient matrix
        mat=matrix(self.A)                                                  # defining co-efficient matrix
        self.b=matrix(self.b).transpose()                                   # defining constant matrix
        L=[[0 for i in range(n)] for j in range(n)]                         # initialising lower triangular matrix
        for i in range(n):
            for j in range(i,n):
                L[j][i]=self.A[j][i]                                        # filling up lower triangular matrix
        L=matrix(L)                                                                   
        U=mat-L                                                             # defining strict upper triangular matrix
        T=((-1)*inv(L))*U                                                   # calculating -U*(inverse of L)     
        C=inv(L)*self.b                                                     # calculating (inverse of L)*b
        diff,count=1,0                                                      # count variable counts number of iterations     
        self.result=matrix([[i] for i in self.initial_guess])                    # initialising solution matrix
        while(diff>accuracy and count<max_iterations):                      # iteration loop
            x_new=T*self.result+C                                           # updating x
            diff=(abs(x_new-self.result)).max()                             # finding accuracy
            self.result=x_new
            count+=1
        if(diff>accuracy):
            return("Iterations exceed limit")                               # maximum iterations exceeded 
        else:
            return(((self.result.transpose()).ravel()).tolist())            # return solution matrix as a list

lss=LinearSystemSolver()                                                    # object creation
for method in ["gauss","gauss-jordan","gauss-siedel"]:
    solution=lss.solve([[4,-1,-1],[-2,6,1],[-1,1,7]],[3,9,-6],method,10000,0.001,[0,0,0])
    print(solution)

# New Class #

class PolynomialSolver:

    # Class PolynomialSolver created by --
    # Name  : SUMAN KUMAR
    # S.no. : 67

    def __init__(self):
        self.order=None
        self.co_eff=None
        self.interval=None
        self.initial_guess=None
        self.max_iterations=None
        self.accuracy=None

    def solve(self,order,co_eff,interval,method,initial_guess1=[-1,1],initial_guess2=1,max_iterations=10000,accuracy=0.0001):

        # initial_guess1 is a list of two starting values for Secant and SecantRF methods, e.g.- [1,3]
        # initial_guess2 is the starting value for NewtonRaphson method, e.g.- 1 
        # initial_guess1, initial_guess2, max_iterations and accuracy are optional arguments

        if(method=="bisection"):
            return (self.BisectionSearch(co_eff,interval,max_iterations,accuracy))
        elif(method=="secant"):
            return (self.Secant(co_eff,initial_guess1,max_iterations,accuracy))
        elif(method=="secantrf"):
            return (self.SecantRF(co_eff,initial_guess1,max_iterations,accuracy))
        else:
            return (self.NewtonRaphson(co_eff,initial_guess2,max_iterations,accuracy))

    def BisectionSearch(self,L,interval,max_iterations,accuracy):           # BisectionSearch function takes 5 arguments
        
        # L is the co-efficient list, e.g.,-
        # If the polynomial is x^2-7x+10, then-
        # L=[10,-7,1]
        # interval is the Interval in which the roots are searched, e.g-
        # interval=[3,6] ; 3-lower limit, 6-upper limit
        # max_iterations is the maximum allowed number of iterations
        # accuracy is the desired Accuracy
        
        self.co_eff=L
        self.interval=interval
        self.max_iterations=max_iterations
        self.accuracy=accuracy
        from numpy.polynomial import polynomial as P
        count=0                                                             # variable to count number of iterations
        while(count<self.max_iterations):
            if(P.polyval(self.interval[0],self.co_eff)*P.polyval(self.interval[1],self.co_eff)<0 and self.interval[1]-self.interval[0]<=self.accuracy):  
                return (self.interval)                                      # return interval (terminal condition)
            m=(self.interval[0]+self.interval[1])/2
            if(P.polyval(self.interval[0],self.co_eff)*P.polyval(m,self.co_eff)<0):
                self.interval=[self.interval[0],m]                          # bisection
            else:
                self.interval=[m,self.interval[1]]                          # bisection
            count+=1
        return ("Iterations exceed limit")
    
    def Secant(self,L,initial_guess,max_iterations,accuracy):               # Secant function takes 5 arguments
        
        # L is the co-efficient list, e.g.,-
        # If the polynomial is x^2-7x+10, then-
        # L=[10,-7,1]
        # initial_guess is a list containing the two starting values, e.g.-
        # initial_guess=[4,6]
        # max_iterations is the maximum allowed number of iterations
        # accuracy is the desired Accuracy
                
        self.co_eff=L
        self.initial_guess=initial_guess
        self.max_iterations=max_iterations
        self.accuracy=accuracy
        from numpy.polynomial import polynomial as P
        x,y=self.initial_guess[0],self.initial_guess[1]                     # x and y store the initial guesses
        count=0                                                             # variable to count number of iterations
        while(count<self.max_iterations):
            z=y-P.polyval(y,self.co_eff)*((y-x)/(P.polyval(y,self.co_eff)-P.polyval(x,self.co_eff))) # calculating next value
            if(abs(P.polyval(z,self.co_eff))<=self.accuracy):               # terminal condition
                return (z)                                                  # return value
            x,y=y,z                                                         # preparing for next iteration
            count+=1
        return ("Iterations exceed limit")

    def SecantRF(self,L,initial_guess,max_iterations,accuracy):             # SecantRF function takes 5 arguments
        
        # L is the co-efficient list, e.g.,-
        # If the polynomial is x^2-7x+10, then-
        # L=[10,-7,1]
        # initial_guess is a list containing the two starting values, e.g.-
        # initial_guess=[1,3]
        # max_iterations is the maximum allowed number of iterations
        # accuracy is the desired Accuracy
        
        self.co_eff=L
        self.initial_guess=initial_guess
        self.max_iterations=max_iterations
        self.accuracy=accuracy
        from numpy.polynomial import polynomial as P
        x,y=self.initial_guess[0],self.initial_guess[1]                     # x and y store the initial guesses
        count=0                                                             # variable to count number of iterations
        while(count<self.max_iterations):
            z=y-P.polyval(y,self.co_eff)*((y-x)/(P.polyval(y,self.co_eff)-P.polyval(x,self.co_eff))) # calculating next value
            if(abs(P.polyval(z,self.co_eff))<=self.accuracy):               # terminal condition
                return (z)                                                  # return value
            if(P.polyval(x,self.co_eff)*P.polyval(z,self.co_eff)>0):        # preparing for next iteration
                x,y=z,y
            else:
                x,y=x,z                                                     # preparing for next iteration
            count+=1
        return ("Iterations exceed limit")

    def NewtonRaphson(self,L,initial_guess,max_iterations,accuracy):        # NewtonRaphson function takes 5 arguments
        
        # L is the co-efficient list, e.g.,-
        # If the polynomial is x^2-9x+14, then-
        # L=[14,-9,1]
        # initial_guess is the starting value, e.g.-
        # initial_guess=8
        # max_iterations is the maximum allowed number of iterations
        # accuracy is the desired Accuracy
        
        self.co_eff=L
        self.initial_guess=initial_guess
        self.max_iterations=max_iterations
        self.accuracy=accuracy
        from numpy.polynomial import polynomial as P
        derivative=P.polyder(self.co_eff)                                   # calculating derivative of polynomial
        x=self.initial_guess                                                # variable to store initial guess
        count=0                                                             # variable to count iterations
        while(count<self.max_iterations):
            y=x-(P.polyval(x,self.co_eff)/P.polyval(x,derivative))          # calculating next value
            if(abs(P.polyval(y,self.co_eff))<=self.accuracy):               # terminal condition
                return (y)                                                  # return value
            x=y                                                             # preparing for next iteration
            count+=1
        return ("Iterations exceed limit")

ps=PolynomialSolver()                                                     # object creation
for method in ["bisection","secant","secantrf","newtonraphson"]:
    solution=ps.solve(2,[14,-9,1],[1,2.5],method,[1,3],1,10000,0.0001)
    print(solution)

# New Class #

class Interpolate:
    
    # Class Interpolate created by --
    # Name  : ROBIN CHAWLA
    # S.no. : 57

    def __init__(self):
        self.x_values=None
        self.fx_values=None
        self.polynomial=None

    def solve(self,L,M,method):
        if(method=="newton"):
            return (self.Newton(L,M))
        else:
            return (self.Lagrange(L,M))

    def Lagrange(self,L,M):                                                 # Lagrange function takes two lists as argument
        
        # L contains the list of x values
        # M contains the list of f(x) values
        # e.g.-
        #L=[1,2,3] , M=[0,-1,0]
        # i.e., f(1)=0, f(2)=-1, f(3)=0
        
        self.x_values=L
        self.fx_values=M
        from numpy import array
        from numpy.polynomial import polynomial as P
        n=len(self.x_values)                                                # n=length of L, i.e., the number of points
        w=(-1*self.x_values[0],1)                                           # initialising polynomial w
        for i in range(1,n):
            w=P.polymul(w,(-1*self.x_values[i],1))                          # calculating w
        result=array([0.0 for i in range(len(w)-1)])                        # initialising result array
        derivative=P.polyder(w)                                             # derivative of w
        for i in range(n):
            result+=(P.polydiv(w,(-1*self.x_values[i],1))[0]*self.fx_values[i])/P.polyval(self.x_values[i],derivative)# calculating result
        result=list(result)                                                 # list of co-efficients
        self.polynomial=""                                                  # string to store final polynomial
        for i in range(len(result)-1,0,-1):                                 # building up the string
            if(result[i]!=0):
                if(result[i]>0 and i!=(len(result)-1)):
                    self.polynomial+=" + "+str(result[i])+"x^"+str(i)+" "
                elif(result[i]>0 and i==(len(result)-1)):
                    self.polynomial+=str(result[i])+"x^"+str(i)+" "
                else:
                    self.polynomial+=" - "+str(-1*result[i])+"x^"+str(i)+" "
        if(result[0]!=0):
            self.polynomial+=" + "+str(result[0]) if result[0]>0 else " - "+str(-1*result[0])
        return (self.polynomial)                                            # return result

    def Newton(self,L,M):                                                   # Newton function takes two lists as arguments

        # L contains the list of x values
        # M contains the list of f(x) values
        # e.g.-
        #L=[1,2,3] , M=[0,-1,0]
        # i.e., f(1)=0, f(2)=-1, f(3)=0
        
        self.x_values=L
        self.fx_values=M
        from numpy import array
        from numpy.polynomial import polynomial as P
        n=len(self.x_values)                                                # n=length of L, i.e., the number of points
        mat=[[0.0 for i in range(n)] for j in range(n)]                     # initialising an n*n matrix 
        for i in range(n):                                                  # filling 1st column of matrix with f(x) values
            mat[i][0]=self.fx_values[i]
        for i in range(1,n):                                                # calculating entries of matrix
            for j in range(n-i):
                mat[j][i]=(mat[j+1][i-1]-mat[j][i-1])/(self.x_values[j+i]-self.x_values[j])
        # The matrix is of the form (for 4 points - x,y,z,w)
        #    f(x)    f(x,y)    f(x,y,z)    f(x,y,z,w)
        #    f(y)    f(y,z)    f(y,z,w)    0
        #    f(z)    f(z,w)    0           0 
        #    f(w)    0         0           0
        
        result=array((mat[0][0],))                                          # initialising result array
        for i in range(1,n):
            prod=(-1*self.x_values[0],1)                                    # initialising prod polynomial which is to be multiplied
                                                                            # with corresponding element of matrix mat
            for j in range(1,i):
                prod=P.polymul(prod,(-1*self.x_values[j],1))                # calculating prod    
            result=P.polyadd(result,array(prod)*mat[0][i])                  # calculating result
        result=list(result)                                                 # list of co-efficients
        self.polynomial=""                                                  # string to store final polynomial
        for i in range(len(result)-1,0,-1):                                 # building up the string
            if(result[i]!=0):
                if(result[i]>0 and i!=(len(result)-1)):
                    self.polynomial+=" + "+str(result[i])+"x^"+str(i)+" "
                elif(result[i]>0 and i==(len(result)-1)):
                    self.polynomial+=str(result[i])+"x^"+str(i)+" "
                else:
                    self.polynomial+=" - "+str(-1*result[i])+"x^"+str(i)+" "
        if(result[0]!=0):
            self.polynomial+=" + "+str(result[0]) if result[0]>0 else " - "+str(-1*result[0])
        return (self.polynomial)                                            # return result

apx=Interpolate()                                                          # object creation
for method in ["newton","lagrange"]:
    solution=apx.solve([-1,0,1,2],[-5,1,3,7],method)
    print(solution)

# New Class #

import math
from math import *
class Integrate:

    # Class PolynomialSolver created by --
    # Name  : D Aravind
    # S.no. : 24

    def __init__(self):
        self.function=None
        self.interval=None
        self.no_of_partitions=None
        self.result=None

    def solve(self,S,interval,method,N=10000):
        
        # N is an optional argument used for controlling the number of partitions
        
        if(method=="trapezoid"):
            return (self.TrapezoidalRule(S,interval,N))
        else:
            return (self.SimpsonsRule(S,interval,N))
    
    def TrapezoidalRule(self,S,interval,N):                                 # The function takes 3 user inputs
        
        # S is the string containing the function to be
        # integrated. The string must contain only those 
        # functions which are available in "math" module,
        # and the syntax used in S must be the same as is
        # used in Python. e.g.,
        # If f(x)=sin(x)^2 + cos(x) + e^(x^0.5) + log(x), then --
        # S=sin(x)**2+cos(x)+exp(x**0.5)+log(x) OR
        # S=pow(sin(x),2)+cos(x)+exp(pow(x,0.5))+log(x)
        # L is a list of two values containing the lower
        # and upper limits of integration. e.g.,-
        # If we have to integrate f(x) from 1.5 to 6, then--
        # L=[1.5,6]
        # N is the variable used for controlling the number
        # of partitions.
        
        self.interval=interval
        self.no_of_partitions=N
        try:
            (lambda x:eval(S))(self.interval[0])
            self.function=lambda x:eval(S)                                  # lambda function which converts user input into a function
        except:
            return("Input is not in proper Python syntax")
        length=(self.interval[1]-self.interval[0])/self.no_of_partitions    # variable to store the length of each sub-interval
        sum_of_values=self.function(self.interval[0])                       # variable to store sum of function values
        for i in range(1,N):                                                # calculating sum of function values
            sum_of_values+=2*self.function(self.interval[0]+length*i)
        sum_of_values+=self.function(self.interval[1])
        self.result=((self.interval[1]-self.interval[0])*sum_of_values)/(2*self.no_of_partitions)  # calculating result
        return (self.result)                                                # return result
    
    def SimpsonsRule(self,S,interval,N):                                    # The function takes 3 user inputs

        # S is the string containing the function to be
        # integrated. The string must contain only those 
        # functions which are available in "math" module,
        # and the syntax used in S must be the same as is
        # used in Python. e.g.,
        # If f(x)=sin(x)^2 + cos(x), then --
        # S=sin(x)**2+cos(x)+exp(x**0.5)+log(x) OR
        # S=pow(sin(x),2)+cos(x)+exp(pow(x,0.5))+log(x)
        # L is a list of two values containing the lower
        # and upper limits of integration. e.g.,-
        # If we have to integrate f(x) from 1.5 to 6, then--
        # L=[1.5,6]
        # N is the variable used for controlling the number
        # of partitions.

        self.interval=interval
        self.no_of_partitions=N
        try:
            (lambda x:eval(S))(self.interval[0])
            self.function=lambda x:eval(S)                                  # lambda function which converts user input into a function
        except:
            return("Input is not in proper Python syntax")
        length=(self.interval[1]-self.interval[0])/(2*self.no_of_partitions)    # variable to store the length of each sub-interval
        sum_of_values=self.function(self.interval[0])                       # variable to store sum of function values
        for i in range(1,2*N):                                              # calculating sum of function values
            temp=self.function(self.interval[0]+length*i)
            sum_of_values+=2*temp if i%2==0 else 4*temp
        sum_of_values+=self.function(self.interval[1])
        self.result=((self.interval[1]-self.interval[0])*sum_of_values)/(6*self.no_of_partitions)      # calculating result
        return (self.result)                                                # return result      
    
igr=Integrate()                                                     # object creation
for method in ["trapezoid","simpson"]:
    solution=igr.solve("sin(x**0.5)+cos(x)+exp(pow(x,2))+log(x)",[1,2],method,10000)
    print (solution)
