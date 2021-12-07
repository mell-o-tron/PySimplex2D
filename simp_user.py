from simpModule import simpView

# Objective:     max c*x
# Constraints:   Ax <= b


A = [[1,-1], [0,1], [1,1], [-1,1], [-1,0]]
B = [4,5];
b = [1, 2, 3, 1, 0]
c = [2,1]

simpView(A, B, b, c)
