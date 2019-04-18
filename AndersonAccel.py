import numpy
import math

def AA_gen(x0, f, m, epsilon=10**-10):
    x1 = f(x0)
    X = [x0, x1]
    F = [x1 - x0]
    k = 1
    mk = 0
    while numpy.linalg.norm(F[-1], 2) > epsilon:
        if k > m:
            mk = m
        else:
            mk = k

        f_k = f(X[-1]) - X[-1]
        F.append(f_k)
        F_mod = [[0. for x in range(mk+1)] for x in range(mk+1)]
        (r, c) = F[-1].shape
        '''
        for i in range(k+1):
            print '[ ',
            for x in range(r):
                print (float)((F[i])[x]),
                print ' ',
            print ' ]'
        '''
        for i in range(mk+1):
            for j in range(mk+1):
                if i == mk and j == mk:
                    F_mod[i][j] = 0.
                elif j == mk:
                    F_mod[i][j] = -.5
                elif i == mk:
                    F_mod[i][j] = 1.
                else:
                    sum = 0.
                    for h in range(r):
                        sum += (float)(F[-j-1][h])*(float)(F[-i-1][h])
                    F_mod[i][j] = sum
        F_mod = numpy.matrix(F_mod,dtype=numpy.float64)
        sol = [0.] * (mk+1)
        sol[mk] = 1.
        sol = numpy.transpose(numpy.matrix(sol,dtype=numpy.float64))

        alphas = numpy.linalg.solve(F_mod, sol)
        x_next = 0.
        (n,r) = alphas.shape
        for i in range(n-1):
            x_next += f(X[-i-1])*(float)(alphas[i])
        X.append(x_next)
        k += 1
        print 'Iteration {}: \n{}'.format(k, X[-1])
    
    print 'Accelerator finished at {} iterations'.format(k)
    return X[-1]
        

'''
Depreciated Anderson Acceleration function with m = 3
Please use AA_gen for a general Anderson Accelerator for
any desired m value
'''
def AA3(x0, f, epsilon=10**-10):
    #Anderson Acceleration for m = 3
    x1 = f(x0)
    X = [x0, x1]
    F = [x1 - x0]
    k = 1
    mk = 0
    alphas = []
    while numpy.linalg.norm((F[-1]),2) > epsilon:
        if k > 3:
            mk = 3
        else:
            mk = k
        #X.append(x_k)
        f_k = f(X[-1]) - X[-1]
        F.append(f_k)
        '''
        print len(F)
        (r, c) = F[0].shape
        for i in range(len(F)):
            print '[ ',
            for x in range(r):
                print (float)(F[i][x]),
                print ' ',
            print ']'
        '''
        a2 = 1.
        a1 = 0.
        a0 = 0.
        if mk == 2:
            #fill in later
            a0 = 0.

            m11 = 0. #x^2
            m12 = 0. #x*y
            m22 = 0. #y^2

            (n, r) = F[-1].shape
            for i in range(n):
                xi = (float)(F[-2][i])
                yi = (float)(F[-1][i])

                m11 += math.pow(xi, 2)
                m22 += math.pow(yi, 2)

                m12 += xi * yi

            mini = numpy.matrix([[m11, m12, -.5],
                                     [m12, m22, -.5],
                                     [1. , 1. , 0  ]],dtype=numpy.float64)
            #print mini
            sol = numpy.matrix([[0.],[0.],[1.]],dtype=numpy.float64)

            alphas = numpy.linalg.solve(mini, sol)
            a1 = (float)(alphas[0])
            a2 = (float)(alphas[1])

        elif mk == 3:
            m11 = 0. #x^2
            m12 = 0. #x*y
            m13 = 0. #x*z
            m22 = 0. #y^2
            m23 = 0. #y*z
            m33 = 0. #z^2
            (n, r) = F[-1].shape
            for i in range(n):
                xi = (float)(F[-3][i])
                yi = (float)(F[-2][i])
                zi = (float)(F[-1][i])

                m11 += math.pow(xi, 2)
                m22 += math.pow(yi, 2)
                m33 += math.pow(zi, 2)

                m12 += xi * yi
                m13 += xi * zi

                m23 += yi * zi

            minimizer=numpy.matrix([[ m11, m12, m13, -.5 ],
                                    [ m12, m22, m23, -.5 ],
                                    [ m13, m23, m33, -.5 ],
                                    [ 1. , 1. , 1. , 0   ]],dtype=numpy.float64)
            #print minimizer
            sol = numpy.transpose(numpy.matrix([0.,0.,0.,1.],dtype=numpy.float64))
            
            alphas = numpy.linalg.solve(minimizer, sol)
            
            a0 = (float)(alphas[0])
            a1 = (float)(alphas[1])
            a2 = (float)(alphas[2])
            #alphas[3] is the lambda value, not needed
        #print alphas
        next_x = 0.
        if len(X) >= 3:
            next_x += a0*f(X[-3])
        if len(X) >= 2:
            next_x += a1*f(X[-2])
        next_x += a2*f(X[-1])
        X.append(next_x)
        print 'Iteration {}: \n{}'.format(k, X[-1])
        k += 1

    print 'AA3 complete after {} iterations'.format(k)
    return X[-1]

