import numpy as np
from sklearn.metrics.pairwise import euclidean_distances

def mixed_integer_euclidean_distances(xp, integer_indices):
    xp = xp.copy()
    xp[:, integer_indices] = np.round(xp[:, integer_indices])
    return euclidean_distances(xp)

def calcRHS(U,d2):
    if U.ndim == 1:
        rhs = np.append(U,np.zeros(d2))
    elif U.ndim == 2:
        rhs = np.vstack((U,np.array([d2*[0]]).T))
    else:
        raise ValueError('U is neither vector nor matirx!')
    
    return rhs

def svdInv(M):

    eps = 1E-14
    u,s,v = np.linalg.svd(M)
    invD = 1/s
    invD[abs(s/s[0])<eps] = 0
    invM = np.matmul(v.T,np.matmul(np.diag(invD),u.T))
    return(invM)


def trainRBF(phi, U, ptail=True, squares=False, xp=None, rho=0.0):
    '''
    Fit cubic RBF interpolation to training data for d>1.
    
    The model for a point z=(z_1,...,z_d) is fitted using n sample points x_1, ..., x_n 
       s(z) = \lambda_1*\Phi(||z-x_1||)+... +\lambda_n*\Phi(||z-x_n||)
                     + c_0 + c_1*z_1 + ... + c_d*z_d  
    
    where \Phi(r)=r^3 denotes the cubic radial basis function. The coefficients \lambda_1, 
    ..., \lambda_n, c_0, c_1, ..., c_d are determined by this training procedure.
    This is for the default case squares==FALSE. In case squares==TRUE 
    there are d additional pure square terms and the model is
    
       s_sq(z) = s(z) + c_d+1*z_1^2 + ... + c_d+d*z_d^2  
    
      
    The linear equation system is solved via SVD inversion. Near-zero elements 
    in the diagonal matrix D are set to zero in D^-1. This is numerically stable 
    for rank-deficient systems.
    
    @param xp      n points x_i of dimension d are arranged in (n x d) matrix xp
    @param U       vector of length n, containing samples u(x_i) of 
                   the scalar function u to be fitted 
                   - or - 
                   (n x m) matrix, where each column 1,...,m contains one vector of samples
                   u_j(x_i) for the m'th model, j=1,...,m
    @param squares [FALSE] flag, see description
    @param rho     [0.0] experimental: 0: interpolating, >0, approximating (spline-like) 
                   Gaussian RBFs
                   
    @return rbf.model,  an object of class RBFinter, which is basically a list 
    with elements:
         coef  (n+d+1 x m) matrix holding in column m the coefficients for the m'th 
                       model:      \lambda_1, ..., \lambda_n, c_0, c_1, ..., c_d.  
                       In case squares==TRUE it is an (n+2d+1 x m) matrix holding  
                       additionally the coefficients c_d+1, ..., c_d+d.                    
         xp  matrix xp   
         d  dimension d 
         npts  number n of points x_i 
         squares  TRUE or FALSE  
         type  "CUBIC"
         
    @seealso   trainGaussRBF, predict.RBFinter, interpRBF
    '''    

    npts = len(xp)
    d2 = None
    pMat = None
    
    if ptail:
        #linear tail
        e = np.array(npts*[1.0])
        pMat = np.column_stack((e,xp))
    if squares:
        # plus direct squares x1**2 x2**2,....
        pMat = np.column_stack((pMat,xp*xp))
    
    d2 = len(pMat[0])
    nMat = np.zeros((d2,d2))
    M = np.column_stack((phi, pMat))
    QQ = pMat.transpose()
    M = np.vstack((M,np.column_stack((QQ,nMat))))
    
    invM = svdInv(M)
    rhs = calcRHS(U,d2)
    coef = np.matmul(invM,rhs)
    
    rbfmodel = dict()
    rbfmodel['coef'] = coef
    rbfmodel['xp'] = xp
    rbfmodel['d'] = d2
    rbfmodel['npts'] = npts
    rbfmodel['ptail'] = ptail
    rbfmodel['squares'] = squares
    rbfmodel['type'] = "CUBIC"     
    return rbfmodel

def trainCubicRBF(xp, U, integer_indices=None, ptail=True, squares=False, rho=0.0):
    if integer_indices is None:
        integer_indices = []
    edist = mixed_integer_euclidean_distances(xp, integer_indices)
    phi = edist*edist*edist
    
    rbfmodel = trainRBF(phi, U, ptail, squares, xp, rho)
    rbfmodel['integer_indices'] = integer_indices
    return rbfmodel

def distLine(x, xp, integer_indices=None):
    z = np.outer(np.ones(len(xp)), x) - xp
    if integer_indices is not None:
        z[:, integer_indices] = np.round(z[:, integer_indices])
    z = np.sqrt(np.sum(z*z, axis=1))
    return z

def interpRBF(x, rbfModel):
    '''
    Apply cubic or Gaussian RBF interpolation to new data for d>1.
    
    param x         vector holding a point of dimension d
    param rbf.model trained RBF model (or set of models), see trainCubicRBF
                     or trainGaussRBF
                   
    return          value s(x) of the trained model at x
                     - or - 
                     vector s_j(x) with values for all trained models j=1,...,m at x
    
    seealso   trainCubicRBF, predict.RBFinter
    '''
    if x.shape[1] != len(rbfModel['xp'][0]):
        raise ValueError('Problem in interpRBF, length of vector and rbf model do not match')
    
    ed = distLine(x, rbfModel['xp'], rbfModel['integer_indices'])  # Euclidean distance of x to all xp, ed is a vector of length nrow(xp)
    
    if rbfModel['type'] == 'CUBIC':
        ph = ed*ed*ed
    else:
        raise ValueError('CUBIC NOT CHOSEN, GAUSS to be implemented')
    
    if rbfModel['ptail']:
        if rbfModel['squares']:
            lhs = np.append(ph, 1)
            lhs = np.append(lhs, x)
            lhs = np.append(lhs, np.multiply(x, x))
        else:
            lhs = np.append(ph, 1)
            lhs = np.append(lhs, x)
    else:
        lhs = ph
    val = np.matmul(lhs, rbfModel['coef'])
    return val

def predictRBFinter(rbfModel, newdata):
    '''
    Apply cubic or Gaussian RBF interpolation
     
    Apply cubic or Gaussian RBF interpolation to a set of new data points for d>1.
    
    param rbf.model trained RBF model (or set of models), see trainCubicRBF 
                     or trainGaussRBF
    param newdata   matrix or data frame with d columns. Each row contains a data point 
                      x_i, i=1,...,n
                    
    return          vector of model responses s(x_i), one element for each data point x_i
                     - or - 
                     if rbf.model is a set of m models, a (n x m)-matrix 
                      containing in each row the response s_j(x_i) of all models 
                     j = 1,...,m to x_i
     
    seealso   trainCubicRBF, trainGaussRBF, interpRBF
    '''

    val = [interpRBF(i, rbfModel) for i in newdata]
    return(val)


def run_rbf_surr():
    xp = np.array([[1, 2.5, 3], [2, 3.5, 4], [3, 4.5, 5]])
    U = np.array([10, 20, 30])
    integer_indices = [0, 2]

    rbfmodel = trainCubicRBF(xp, U, integer_indices)
    print(rbfmodel)
    newdata = np.asmatrix([[3, 3.5, 4]])
    predictions = predictRBFinter(rbfmodel, newdata)

    print(predictions)
