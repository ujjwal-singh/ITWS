class LPsolver:

    # Class LPsolver created by --
    # Name  : UJJWAL SINGH
    # S.no. : 70

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
        
        no_of_variables=len(L)                                              # variable to store total number of variables
        no_of_equations=len(M)                                              # variable to store total number of equations(constraints)
        mat=[[0 for i in range(no_of_variables+no_of_equations+2)] for j in range(no_of_equations+1)]   # initialising the augmented matrix
        for i in range(no_of_equations):
            for j in range(no_of_variables):
                mat[i][j]=M[i][j]                                           # filling up the matrix
        for i in range(no_of_variables):
            mat[no_of_equations][i]=(-1)*L[i]                               # filling up the matrix
        for i in range(no_of_equations+1):
            mat[i][no_of_variables+i]=1                                     # filling up the matrix
        for i in range(no_of_equations):
            mat[i][no_of_variables+no_of_equations+1]=M[i][-1]              # filling up the matrix
        mat[-1][-1]=0                                                       # matrix ready
        for i in range(no_of_equations+1):
            mat[i][no_of_variables+i]=mat[i][-1]/mat[i][no_of_variables+i]
        flag=0
        if(min(mat[no_of_equations])>=0):                                   # testing terminal condition
            flag=1
        while(flag==0):                                                     # loop
            pivot_col=mat[no_of_equations].index(min(mat[no_of_equations])) # finding column of pivot
            pivot_row,temp,flag2=-1,0,0                                         
            for i in range(no_of_equations):                                                        
                if(mat[i][pivot_col]>0 and (pivot_row==-1 or mat[i][-1]/mat[i][pivot_col]<temp)):   # finding row of pivot
                    pivot_row=i
                    temp=mat[i][-1]/mat[i][pivot_col]
                    flag2=1
            if(flag2==0):
                return ("Unbound")                                          # terminal condition 
            pivot=mat[pivot_row][pivot_col]
            for i in range(no_of_equations+1):                              # loop to set all other elements of pivot's column to 0 by row operations
                if(i==pivot_row):
                    continue
                temp=mat[i][pivot_col]/pivot
                for j in range(no_of_variables+no_of_equations+2):
                    if(j==pivot_col):
                        mat[i][j]=0.0
                    else:    
                        mat[i][j]=mat[i][j]-(temp*mat[pivot_row][j])
            if(min(mat[no_of_equations])>=0):
                flag=1                                                      # terminal condition
        st=""                                                               # string to store final answer
        for i in range(no_of_variables):                                    # building up the string
            count,temp,flag=0,0,1
            for j in range(no_of_equations+1):
                if(mat[j][i]!=0):
                    count+=1
                    if(count>1):
                        flag=0
                        break
                    temp=mat[j][-1]/mat[j][i]
            if(flag==1):
                st+="x_"+str(i+1)+" = "+str(temp)+" ; "
            else:
                st+="x_"+str(i+1)+"=0 ; "
        st+="Maximum value of objective = "+str(mat[-1][-1])
        return (st)                                                         # return string

lps=LPsolver()                                                              # object creation
solution=lps.solve([1,3],[[1,2,10],[1,0,5],[0,1,4]],option="Simplex")     # calling Simplex method
print(solution)

# New Class #

class LinearSystemSolver:

    # Class LinearSystemSolver created by --
    # Name  : SUSHISH KUMAR
    # S.no. : 69

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
        
        n=len(A)                                                            # finding order of co-efficient matrix
        mat=[[0 for i in range(n+1)] for j in range(n)]                     # initialising the augmented matrix
        for i in range(n):
            for j in range(n):
                mat[i][j]=A[i][j]                                           # filling up the entries of matrix A
        for i in range(n):
            mat[i][n]=B[i]                                                  # filling up the entries of matrix B
        for i in range(n-1):
            if(mat[i][i]==0):                                               # checking if the pivot is 0
                flag=0  
                for j in range(i+1,n):
                    if(mat[j][i]!=0):
                        temp=mat[j]                                         # swapping rows if pivot is 0
                        mat[j]=mat[i]
                        mat[i]=temp
                        flag=1
                        break
                if(flag==0):                                                
                    continue                                                # continue if the column contains only 0s
            for j in range(i+1,n):
                temp=mat[j][i]/mat[i][i]                                    
                for k in range(n+1):                                        # performing row operations to make the entries below the pivot=0
                    if(k==i):
                        mat[j][k]=0.0
                    else:
                        mat[j][k]=mat[j][k]-temp*mat[i][k]
        solutions=[0 for i in range(n)]                                     # initialisng the matrix containing the solution matrix
        for i in range(n-1,-1,-1):
            s=0
            for j in range(n):                                              # solving using back substitution
                s+=mat[i][j]*solutions[j]
            solutions[i]=(mat[i][n]-s)/mat[i][i]
        return(solutions)                                                   # returning solution matrix

    def GaussJordan(self,A,B):                                              # GaussJordan function takes two lists as arguments
                                                                            # A is the co-efficient matrix, B is matrix of constants

        # Suppose we have to solve the system --
        # 1.x+1.y+1.z=6
        # 2.x+1.y+(-1).z=1
        # (-1).x+2.y+2.z=9
        # Then --
        # A=[[1,1,1],[2,1,-1],[-1,2,2]]
        # B=[6,1,9]
        
        n=len(A)                                                            # finding order of co-efficient matrix
        mat=[[0 for i in range(n+1)] for j in range(n)]                     # initialising the augmented matrix
        for i in range(n):
            for j in range(n):
                mat[i][j]=A[i][j]                                           # filling up the entries of matrix A
        for i in range(n):
            mat[i][n]=B[i]                                                  # filling up the entries of matrix B
        for i in range(n):
            if(mat[i][i]==0):                                               # checking if the pivot is 0
                flag=0
                for j in range(n):
                    if(mat[j][i]!=0):
                        temp=mat[j]                                         # swapping rows if pivot is 0
                        mat[j]=mat[i]
                        mat[i]=temp
                        flag=1
                        break
                if(flag==0):
                    continue                                                # continue if the column contains only 0s
            for j in range(n):
                if(j==i):
                    continue
                temp=mat[j][i]/mat[i][i]
                for k in range(n+1):                                        # performing row operations to make the entries in the the pivot column=0
                    if(k==i):
                        mat[j][k]=0.0
                    else:
                        mat[j][k]=mat[j][k]-temp*mat[i][k]
        solutions=[0 for i in range(n)]                                     # initialisng the matrix containing the solution matrix
        for i in range(n):                                                  
            solutions[i]=mat[i][n]/mat[i][i]                                # solving
        return(solutions)                                                   # returning solution matrix
    
    def GaussSiedel(self,A,B,max_iterations,accuracy,initial_guess):        # GaussSiedel function takes two lists as arguments 
                                                                            # A is the co-efficient matrix, B is matrix of constants         

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
        
        from numpy import matrix                                            # maximmum allowed iterations, desired accuracy and initial guess as arguments
        from numpy.linalg import inv
        n=len(A)                                                            # finding order of co-efficient matrix
        mat=matrix(A)                                                       # defining co-efficient matrix
        b=matrix(B).transpose()                                             # defining constant matrix
        L=[[0 for i in range(n)] for j in range(n)]                         # initialising lower triangular matrix
        for i in range(n):
            for j in range(i,n):
                L[j][i]=A[j][i]                                             # filling up lower triangular matrix
        L=matrix(L)                                                                   
        U=mat-L                                                             # defining strict upper triangular matrix
        T=((-1)*inv(L))*U                                                   # calculating -U*(inverse of L)     
        C=inv(L)*b                                                          # calculating (inverse of L)*b
        diff,count=1,0                                                      # count variable counts number of iterations     
        x=matrix([[i] for i in initial_guess])                              # initialising solution matrix
        while(diff>accuracy and count<max_iterations):                      # iteration loop
            x_new=T*x+C                                                     # updating x
            diff=(abs(x_new-x)).max()                                       # finding accuracy
            x=x_new
            count+=1
        if(diff>accuracy):
            return("Iterations exceed limit")                               # maximum iterations exceeded 
        else:
            return(x.transpose())                                           # return solution matrix

lss=LinearSystemSolver()                                                   # object creation
for method in ["gauss","gauss-jordan","gauss-siedel"]:
    solution=lss.solve([[1,1,1],[2,1,-1],[-1,2,2]],[6,1,9],method,10000,0.001,[1,1,1])
    print(solution)
