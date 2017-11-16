import numpy as np
import matplotlib.pyplot as plt


def bvector(x):
    g = np.array([[0],[0]])
    g = np.transpose(g)
    return(g)


def rk3(A, bvector, y0, interval, N):
    alpha1 = 1000
    alpha2 = 1
    x = np.linspace(interval[0], interval[1], N+1, dtype = 'float')
    n = len(y0)
    h = (interval[1] - interval[0]) / N
    y = np.zeros((n, N+1))
    y[:,0] = np.transpose(y0)

    for i in range(0, N):
        y1 = y[:,i] + h * (np.matmul(A, y[:,i]) + bvector(x[i]))
        print(h * (np.matmul(A, y[:,i]) + bvector(x[i])))
        print(y[:,i])
        y2 = (3/4) * y[:,i] + (1/4) * y1 + (h/4) * (np.matmul(A, y1) + bvector(x[i]+h))
        y[:,i+1] = (1/3) * y[:,i] + (2/3) * y2 + h * (2/3) * (np.matmul(A, y2) + bvector(x[i]+h))
    print(y)
    y_exact = np.zeros_like(y)
    for l in range(len(x)):
        y_exact[0,l] = np.exp(-alpha1 * x[l])
        y_exact[1,l] = (alpha1 / (alpha1 - alpha2)) * (np.exp(-alpha2 * x[l]) - np.exp(-alpha1 * x[l]))
    print(y_exact)
    return(x, y)

rk3(np.array([[-1000, 0],[1000, -1]]), bvector, np.array([[1], [0]]), [0,0.1], 400)

def dirk3(A, bvector, y0, interval, N):
    alpha1 = 1000
    alpha2 = 1
    x = np.linspace(interval[0], interval[1], N+1, dtype = 'float')
    n = len(y0)
    h = (interval[1] - interval[0]) / N
    y = np.zeros((n, N+1))
    y[:,0] = np.transpose(y0)
    
    mu = 0.5 * (1 - (1/np.sqrt(3)))
    nu = 0.5 * (np.sqrt(3) - 1)
    gam = 3 / (2 * (3 + np.sqrt(3)))
    lam = (3 * (1 + np.sqrt(3))) / (2 * (3 + np.sqrt(3)))
    
    for j in range(0, N):
        y1 = np.matmul(y[:,j], np.linalg.inv(np.identity(2) - (h * mu * A))) +\
        np.matmul((h * mu * bvector(x[j] + h * mu)), np.linalg.inv(np.identity(2) - (h * mu * A)))
        y2 = np.matmul(y1, np.linalg.inv(np.identity(2) - (h * mu * A))) + np.matmul()
        print(y1)
        
dirk3(np.array([[-1000, 0],[1000, -1]]), bvector, np.array([[1], [0]]), [0,0.1], 400)
        