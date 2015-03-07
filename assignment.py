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
    
    def solve(self,A,B,method,max_iterations,accuracy,initial_guess):
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
                    if(self.augmented[j][i]!=0):
                        temp=self.augmented[j]                              # swapping rows if pivot is 0
                        self.augmented[j]=self.augmented[i]
                        self.augmented[i]=temp
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
                    if(self.augmented[j][i]!=0):
                        temp=self.augmented[j]                              # swapping rows if pivot is 0
                        self.augmented[j]=self.augmented[i]
                        self.augmented[i]=temp
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
                                                                            # maximmum allowed iterations, desired accuracy and initial guess as arguments

        # Suppose we have to solve the system --
        # 1.x+1.y+1.z=10
        # 0.x+1.y+0.z=2
        # 0.x+0.y+1.z=5
        # Then --
        # A=[[1,1,1],[0,1,0],[0,0,1]]
        # B=[10,2,5]
        # max_iterations is the maximum number of iterations
        # accuracy is the desired amount of Accuracy
        # initial_guess is the starting matrix, e.g.--
        #initial_guess=[1,1,1]
        
        from numpy import matrix                                            
        from numpy.linalg import inv
        self.A=A
        self.b=B
        self.max_iterations=max_iterations
        self.accuracy=accuracy
        self.initial_guess=initial_guess
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
        self.result=matrix([[i] for i in initial_guess])                    # initialising solution matrix
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
    solution=lss.solve([[16,3],[7,-11]],[11,13],method,10000,0.001,[1,1])
    print(solution)

# New Class #

class PolynomialSolver:

    # Class PolynomialSolver created by --
    # Name  : SUMAN KUMAR
    # S.no. : 67

    def solve(self,order,co_eff,method,interval,initial_guess1,initial_guess2,max_iterations,accuracy):

        # initial_guess1 is a list of two starting values for Secant and SecantRF methods, e.g.- [1,3]
        # initial_guess2 is the starting value for NewtonRaphson method, e.g.- 1 

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
        
        from numpy.polynomial import polynomial as P
        count=0                                                             # variable to count number of iterations
        while(count<max_iterations):
            if(P.polyval(interval[0],L)*P.polyval(interval[1],L)<0 and interval[1]-interval[0]<=accuracy):  # terminal condition
                return (interval)                                           # return interval
            m=(interval[0]+interval[1])/2
            if(P.polyval(interval[0],L)*P.polyval(m,L)<0):
                interval=[interval[0],m]                                    # bisection
            else:
                interval=[m,interval[1]]                                    # bisection
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
                
        from numpy.polynomial import polynomial as P
        x,y=initial_guess[0],initial_guess[1]                               # x and y store the initial guesses
        count=0                                                             # variable to count number of iterations
        while(count<max_iterations):
            z=y-P.polyval(y,L)*((y-x)/(P.polyval(y,L)-P.polyval(x,L)))      # calculating next value
            if(abs(P.polyval(z,L))<=accuracy):                              # terminal condition
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
        
        from numpy.polynomial import polynomial as P
        x,y=initial_guess[0],initial_guess[1]                               # x and y store the initial guesses
        count=0                                                             # variable to count number of iterations
        while(count<max_iterations):
            z=y-P.polyval(y,L)*((y-x)/(P.polyval(y,L)-P.polyval(x,L)))      # calculating next value
            if(abs(P.polyval(z,L))<=accuracy):                              # terminal condition
                return (z)                                                  # return value
            if(P.polyval(x,L)*P.polyval(z,L)>0):                            # preparing for next iteration
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
        
        from numpy.polynomial import polynomial as P
        derivative=P.polyder(L)                                             # calculating derivative of polynomial
        x=initial_guess                                                     # variable to store initial guess
        count=0                                                             # variable to count iterations
        while(count<max_iterations):
            y=x-(P.polyval(x,L)/P.polyval(x,derivative))                    # calculating next value
            if(abs(P.polyval(y,L))<=accuracy):                              # terminal condition
                return (y)                                                  # return value
            x=y                                                             # preparing for next iteration
            count+=1
        return ("Iterations exceed limit")

ps=PolynomialSolver()                                                     # object creation
for method in ["bisection","secant","secantrf","newtonraphson"]:
    solution=ps.solve(2,[14,-9,1],method,[1,2.5],[1,3],1,1000,0.00001)
    print(solution)

# New Class #

class Interpolate:
    
    # Class Interpolate created by --
    # Name  : ROBIN CHAWLA
    # S.no. : 57

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
        
        from numpy import array
        from numpy.polynomial import polynomial as P
        n=len(L)                                                            # n=length of L, i.e., the number of points
        w=(-1*L[0],1)                                                       # initialising polynomial w
        for i in range(1,n):
            w=P.polymul(w,(-1*L[i],1))                                      # calculating w
        result=array([0.0 for i in range(len(w)-1)])                        # initialising result array
        derivative=P.polyder(w)                                             # derivative of w
        for i in range(n):
            result+=(P.polydiv(w,(-1*L[i],1))[0]*M[i])/P.polyval(L[i],derivative)   # calculating result
        co_eff=list(result)                                                 # list of co-efficients
        st=""                                                               # string to store final polynomial
        for i in range(len(co_eff)-1,0,-1):                                 # building up the string
            if(co_eff[i]!=0):
                if(co_eff[i]>0 and i!=(len(co_eff)-1)):
                    st+=" + "+str(co_eff[i])+"x^"+str(i)+" "
                elif(co_eff[i]>0 and i==(len(co_eff)-1)):
                    st+=str(co_eff[i])+"x^"+str(i)+" "
                else:
                    st+=" - "+str(-1*co_eff[i])+"x^"+str(i)+" "
        if(co_eff[0]!=0):
            st+=" + "+str(co_eff[0]) if co_eff[0]>0 else " - "+str(-1*co_eff[0])
        return (st)                                                         # return answer

    def Newton(self,L,M):                                                   # Newton function takes two lists as arguments

        # L contains the list of x values
        # M contains the list of f(x) values
        # e.g.-
        #L=[1,2,3] , M=[0,-1,0]
        # i.e., f(1)=0, f(2)=-1, f(3)=0
        
        from numpy import array
        from numpy.polynomial import polynomial as P
        n=len(L)                                                            # n=length of L, i.e., the number of points
        mat=[[0.0 for i in range(n)] for j in range(n)]                     # initialising an n*n matrix 
        for i in range(n):                                                  # filling 1st column of matrix with f(x) values
            mat[i][0]=M[i]
        for i in range(1,n):                                                # calculating entries of matrix
            for j in range(n-i):
                mat[j][i]=(mat[j+1][i-1]-mat[j][i-1])/(L[j+i]-L[j])
        # The matrix is of the form (for 4 points - x,y,z,w)
        #    f(x)    f(x,y)    f(x,y,z)    f(x,y,z,w)
        #    f(y)    f(y,z)    f(y,z,w)    0
        #    f(z)    f(z,w)    0           0 
        #    f(w)    0         0           0
        
        result=array((mat[0][0],))                                          # initialising result array
        for i in range(1,n):
            prod=(-1*L[0],1)                                                # initialising prod polynomial which is to be multiplied
                                                                            # with corresponding element of matrix mat
            for j in range(1,i):
                prod=P.polymul(prod,(-1*L[j],1))                            # calculating prod    
            result=P.polyadd(result,array(prod)*mat[0][i])                  # calculating result
        co_eff=list(result)                                                 # list of co-efficients
        st=""                                                               # string to store final polynomial
        for i in range(len(co_eff)-1,0,-1):                                 # building up the string
            if(co_eff[i]!=0):
                if(co_eff[i]>0 and i!=(len(co_eff)-1)):
                    st+=" + "+str(co_eff[i])+"x^"+str(i)+" "
                elif(co_eff[i]>0 and i==(len(co_eff)-1)):
                    st+=str(co_eff[i])+"x^"+str(i)+" "
                else:
                    st+=" - "+str(-1*co_eff[i])+"x^"+str(i)+" "
        if(co_eff[0]!=0):
            st+=" + "+str(co_eff[0]) if co_eff[0]>0 else " - "+str(-1*co_eff[0])
        return (st)                                                         # return answer

apx=Interpolate()                                                          # object creation
for method in ["newton","lagrange"]:
    solution=apx.solve([1,2,3],[3,7,13],method)
    print(solution)

# New Class #

from math import *
class Integrate:

    # Class PolynomialSolver created by --
    # Name  : D Aravind
    # S.no. : 24

    def solve(self,S,interval,N,method):
        if(method=="trapezoid"):
            return (self.TrapezoidalRule(S,interval,N))
        else:
            return (self.SimpsonsRule(S,interval,N))
    
    def TrapezoidalRule(self,S,interval,N):                             # The function takes 3 user inputs
        
        # S is the string containing the function to be
        # integrated. The string must contain only those 
        # functions which are available in "math" module,
        # and the syntax used in S must be the same as is
        # used in Python. e.g.,
        # If f(x)=sin(x)^2 + cos(x), then --
        # S=sin(x)**2+cos(x) OR
        # S=pow(sin(x),2)+cos(x)
        # L is a list of two values containing the lower
        # and upper limits of integration. e.g.,-
        # If we have to integrate f(x) from 1.5 to 6, then--
        # L=[1.5,6]
        # N is the variable used for controlling the number
        # of partitions.
        
        fun=lambda x:eval(S)                                            # lambda function which converts user input into a function
        length=(interval[1]-interval[0])/N                              # variable to store the length of each sub-interval
        sum_of_values=fun(interval[0])                                  # variable to store sum of function values
        for i in range(1,N):                                            # calculating sum of function values
            sum_of_values+=2*fun(interval[0]+length*i)
        sum_of_values+=fun(interval[1])
        result=((interval[1]-interval[0])*sum_of_values)/(2*N)          # calculating result
        return (result)                                                 # return result
    
    def SimpsonsRule(self,S,interval,N):                                # The function takes 3 user inputs

        # S is the string containing the function to be
        # integrated. The string must contain only those 
        # functions which are available in "math" module,
        # and the syntax used in S must be the same as is
        # used in Python. e.g.,
        # If f(x)=sin(x)^2 + cos(x), then --
        # S=sin(x)**2+cos(x) OR
        # S=pow(sin(x),2)+cos(x)
        # L is a list of two values containing the lower
        # and upper limits of integration. e.g.,-
        # If we have to integrate f(x) from 1.5 to 6, then--
        # L=[1.5,6]
        # N is the variable used for controlling the number
        # of partitions.

        fun=lambda x:eval(S)                                        # lambda function which converts user input into a function
        length=(interval[1]-interval[0])/(2*N)                      # variable to store the length of each sub-interval
        sum_of_values=fun(interval[0])                              # variable to store sum of function values
        for i in range(1,2*N):                                      # calculating sum of function values
            temp=fun(interval[0]+length*i)
            sum_of_values+=2*temp if i%2==0 else 4*temp
        sum_of_values+=fun(interval[1])
        result=((interval[1]-interval[0])*sum_of_values)/(6*N)      # calculating result
        return (result)                                             # return result      
    
igr=Integrate()                                                     # object creation
for method in ["trapezoid","simpson"]:
    solution=igr.solve("pow(sin(x),2)+cos(x)+x",[0,1.57],10000,method)
    print (solution)
