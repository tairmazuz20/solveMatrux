"""
רעות אביטן
שילת חכימי
תאיר מזוז
"""
def getMatrixMinor(m, i, j):
    return [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]

#Determination of determination#
def getMatrixDet(m):
    # base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1) ** c) * m[0][c] * getMatrixDet(getMatrixMinor(m, 0, c))
    return determinant


#Arrangement of a matrix#
def orderMatrix(a):
    for i in range(len(a)):
        for j in range(len(a)):
            if i == j and a[i][j] == 0:
                q = a [i]
                for n in range(len(a)):
                    if a[n][i] != 0:
                        a[i] = a[n]
                        a[n]=q


#Matrix multiplication 4X4#
def multiMatrix(m1, m2):

    res = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    for i in range(len(m1)):
        # iterate through columns of Y
        for j in range(len(m2[0])):
            # iterate through rows of Y
            for k in range(len(m2)):
                res[i][j] += m1[i][k] * m2[k][j]

    return res


#Matrix multiplication 4X1#
def multiM(m1, m2):

    res1 = [[0,], [0,], [0,], [0,]]
    for i in range(len(m1)):
        # iterate through columns of Y
        for j in range(len(m2[0])):
            # iterate through rows of Y
            for k in range(len(m1)):
                res1[i][j] += (m1[i][k] * m2[k][j])

    return res1

#Inversion matrix calculation#
def reverseMatrix(m1):
    if getMatrixDet(m1)!=0:
        orderMatrix(m1)
        IA = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
        res = multiMatrix(IA, IA)
        for j in range(len(m1)):
            for i in range(len(m1)):
                if i >= j:
                    IA = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
                    if i == j:
                        dev = m1[i][j]
                        for k in range(len(m1)):
                            m1[i][k] /= dev
                        IA[i][j] /=dev
                        res = multiMatrix(IA, res)
                    else:
                        if m1[i][j] != 0:
                            mu = -m1[i][j]
                            for k in range(len(m1)):
                                m1[i][k] += mu*m1[j][k]
                            IA[i][j] = mu
                            res = multiMatrix(IA, res)
        for j in range(len(m1)-1,-1,-1):
            for i in range(len(m1)-2,-1,-1):
                if j > i:
                    IA = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
                    if m1[i][j] != 0:
                        mu = -m1[i][j]
                        for k in range(len(m1)):
                            m1[i][k] += mu * m1[j][k]
                        IA[i][j] = mu
                        res = multiMatrix(IA, res)
        return res

#Decomposition into L and U matrices#
def LU():
    m = [[1, 0, 1, 0],
         [1, 5, 0, 0],
         [0, 0, 1, 2],
         [1, 1, 0, 7]]
    IA = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    ll = multiMatrix(IA, IA)
    for j in range(len(m)):
        for i in range(len(m)):
            if i > j:
                IA = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
                if m[i][j] != 0:
                    IA[i][j] = -m[i][j]/m[j][j]
                    m = multiMatrix(IA, m)
                    IA[i][j] = -IA[i][j]
                    ll = multiMatrix(ll,IA)
    return ll,m


#Matrix solution using the LU method#
def solveLU(b):
    ll, uu =LU()
    l1 = reverseMatrix(ll)
    u1 = reverseMatrix(uu)
    y = multiM(l1,b)
    x = multiM(u1, y)
    print("The solve of matrix is:")
    for a in x:
        print(a)

#The main plan#
def main():

    A = [[1, 0, 1, 0],
         [1, 5, 0, 0],
        [0, 0, 1, 2],
        [1, 1, 0, 7]]
    B= [[1],
        [2],
        [3],
        [1]]
    CA = [[1, 0, 1, 0],
         [1, 5, 0, 0],
         [0, 0, 1, 2],
         [1, 1, 0, 7]]
    print("CA")
    for l in CA:
        print(l)
    A1 = reverseMatrix(A)
    print("A1- The reverse Matrix is:")
    for l in A1:
        print(l)
    test = multiMatrix(A1, CA)
    print("test:")
    for l in test:
        print(l)
    solveLU(B)

main()