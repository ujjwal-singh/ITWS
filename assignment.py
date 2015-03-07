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
