# pylmm is a python-based linear mixed-model solver with applications to GWAS
# Copyright (C) 2015  Nicholas A. Furlotte (nick.furlotte@gmail.com)

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np, matplotlib.pyplot as plt, scipy.stats as ss
from scipy import linalg
from scipy import optimize
from scipy import stats


def matrixMult(A, B):
    # If there is no fblas then we will revert to np.dot()
    try:
        linalg.fblas
    except AttributeError:
        return np.dot(A, B)

    # If the matrices are in Fortran order then the computations will be faster
    # when using dgemm.  Otherwise, the function will copy the matrix and that takes time.
    if not A.flags['F_CONTIGUOUS']:
        AA = A.T
        transA = True
    else:
        AA = A
        transA = False

    if not B.flags['F_CONTIGUOUS']:
        BB = B.T
        transB = True
    else:
        BB = B
        transB = False

    return linalg.fblas.dgemm(alpha=1., a=AA, b=BB, trans_a=transA, trans_b=transB)


def calculateKinship(W, center=False):
    """
    W is an n x m matrix encoding SNP minor alleles.
    This function takes a matrix oF SNPs, imputes missing values with the maf,
    normalizes the resulting vectors and returns the RRM matrix.
    """
    n = W.shape[0]
    m = W.shape[1]
    keep = []
    # for i in range(m):
    #   mn = W[True - np.isnan(W[:,i]),i].mean()
    #   W[np.isnan(W[:,i]),i] = mn
    #   vr = W[:,i].var()
    #   if vr == 0: continue

    #   keep.append(i)
    #   W[:,i] = (W[:,i] - mn) / np.sqrt(vr)

    K = matrixMult(W, W.T) * 1.0 / float(m)
    if center:
        P = np.diag(np.repeat(1, n)) - 1 / float(n) * np.ones((n, n))
        S = np.trace(matrixMult(matrixMult(P, K), P))
        K_n = (n - 1) * K / S
        return K_n
    return K


class LMM:
    """
    This is a simple version of EMMA/fastLMM.
    The main purpose of this module is to take a phenotype vector (Y), a set of covariates (X) and a kinship matrix (K)
    and to optimize this model by finding the maximum-likelihood estimates for the model parameters.
    There are three model parameters: heritability (h), covariate coefficients (beta) and the total
    phenotypic variance (sigma).
    Heritability as defined here is the proportion of the total variance (sigma) that is attributed to
    the kinship matrix.

    For simplicity, we assume that everything being input is a numpy array.
    If this is not the case, the module may throw an error as conversion from list to numpy array
    is not done consistently.
   """

    def __init__(self, Y, K, Kva=[], Kve=[], X0=None, verbose=False):

        """
      The constructor takes a phenotype vector or array of size n.
      It takes a kinship matrix of size n x n.  Kva and Kve can be computed as Kva,Kve = linalg.eigh(K) and cached.
      If they are not provided, the constructor will calculate them.
      X0 is an optional covariate matrix of size n x q, where there are q covariates.
      When this parameter is not provided, the constructor will set X0 to an n x 1 matrix of all ones to represent a mean effect.
      """

        if X0 == None:
            X0 = np.ones(len(Y)).reshape(len(Y), 1)
        self.verbose = verbose

        x = True - np.isnan(Y)
        x = x.reshape(-1, )
        if not x.sum() == len(Y):
            if self.verbose: print("Removing %d missing values from Y\n" % ((True - x).sum()))
            Y = Y[x]
            K = K[x, :][:, x]
            X0 = X0[x, :]
            Kva = []
            Kve = []
        self.nonmissing = x

        if len(Kva) == 0 or len(Kve) == 0:
            if self.verbose: print(
                "Obtaining eigendecomposition for %dx%d matrix\n" % (K.shape[0], K.shape[1]))
            Kva, Kve = linalg.eigh(K)
            if self.verbose: print("Total time: %0.3f\n" % (end - begin))

        self.K = K
        self.Kva = Kva
        self.Kve = Kve
        self.N = self.K.shape[0]
        self.Y = Y.reshape((self.N, 1))
        self.X0 = X0

        if sum(self.Kva < 1e-6):
            if self.verbose: print("Cleaning %d eigen values\n" % (sum(self.Kva < 0)))
            self.Kva[self.Kva < 1e-6] = 1e-6

        self.transform()

    def transform(self):

        """
        Computes a transformation on the phenotype vector and the covariate matrix.
        The transformation is obtained by left multiplying each parameter by the transpose of the
        eigenvector matrix of K (the kinship).
        """

        self.Yt = matrixMult(self.Kve.T, self.Y)
        self.X0t = matrixMult(self.Kve.T, self.X0)
        self.X0t_stack = np.hstack([self.X0t, np.ones((self.N, 1))])
        self.q = self.X0t.shape[1]

    def getMLSoln(self, h, X):

        """
        Obtains the maximum-likelihood estimates for the covariate coefficients (beta),
        the total variance of the trait (sigma) and also passes intermediates that can
        be utilized in other functions. The input parameter h is a value between 0 and 1 and represents
        the heritability or the proportion of the total variance attributed to genetics.  The X is the
        covariate matrix.
        """

        S = 1.0 / (h * self.Kva + (1.0 - h))
        Xt = X.T * S
        XX = matrixMult(Xt, X)
        XX_i = linalg.pinv(XX)
        beta = matrixMult(matrixMult(XX_i, Xt), self.Yt)
        Yt = self.Yt - matrixMult(X, beta)
        Q = np.dot(Yt.T * S, Yt)
        sigma = Q * 1.0 / (float(self.N) - float(X.shape[1]))
        return beta, sigma, Q, XX_i, XX

    def LL_brent(self, h, X=None, REML=False):
        # brent will not be bounded by the specified bracket.
        # I return a large number if we encounter h < 0 to avoid errors in LL computation during the search.
        if h < 0: return 1e6
        return -self.LL(h, X, stack=False, REML=REML)[0]

    def LL(self, h, X=None, stack=True, REML=False):

        """
        Computes the log-likelihood for a given heritability (h).  If X==None, then the
        default X0t will be used.  If X is set and stack=True, then X0t will be matrix concatenated with
        the input X.  If stack is false, then X is used in place of X0t in the LL calculation.
        REML is computed by adding additional terms to the standard LL and can be computed by setting REML=True.
        """

        if X is None:
            X = self.X0t
        elif stack:
            self.X0t_stack[:, (self.q)] = matrixMult(self.Kve.T, X)[:, 0]
            X = self.X0t_stack

        n = float(self.N)
        q = float(X.shape[1])
        beta, sigma, Q, XX_i, XX = self.getMLSoln(h, X)
        LL = n * np.log(2 * np.pi) + np.log(h * self.Kva + (1.0 - h)).sum() + n + n * np.log(1.0 / n * Q)
        LL = -0.5 * LL

        if REML:
            LL_REML_part = q * np.log(2.0 * np.pi * sigma) + np.log(linalg.det(matrixMult(X.T, X))) - np.log(
                linalg.det(XX))
            LL = LL + 0.5 * LL_REML_part

        LL = LL.sum()
        return LL, beta, sigma, XX_i

    def getMax(self, H, X=None, REML=False):

        """
        Helper functions for .fit(...).
        This function takes a set of LLs computed over a grid and finds possible regions
        containing a maximum.  Within these regions, a Brent search is performed to find the
        optimum.
        """

        n = len(self.LLs)
        HOpt = []
        for i in range(1, n - 2):
            if self.LLs[i - 1] < self.LLs[i] and self.LLs[i] > self.LLs[i + 1]:
                HOpt.append(optimize.brent(self.LL_brent, args=(X, REML), brack=(H[i - 1], H[i + 1])))
                if np.isnan(HOpt[-1]): HOpt[-1] = H[i - 1]
                # if np.isnan(HOpt[-1]): HOpt[-1] = self.LLs[i-1]
                # if np.isnan(HOpt[-1][0]): HOpt[-1][0] = [self.LLs[i-1]]

        if len(HOpt) > 1:
            if self.verbose: print("NOTE: Found multiple optima.  Returning first...\n")
            return HOpt[0]
        elif len(HOpt) == 1:
            return HOpt[0]
        elif self.LLs[0] > self.LLs[n - 1]:
            return H[0]
        else:
            return H[n - 1]

    def fit(self, X=None, ngrids=100, REML=True):

        """
        Finds the maximum-likelihood solution for the heritability (h) given the current parameters.
        X can be passed and will transformed and concatenated to X0t.  Otherwise, X0t is used as
        the covariate matrix.

        This function calculates the LLs over a grid and then uses .getMax(...) to find the optimum.
        Given this optimum, the function computes the LL and associated ML solutions.
        """

        if X is None:
            X = self.X0t
        else:
            # X = np.hstack([self.X0t,matrixMult(self.Kve.T, X)])
            self.X0t_stack[:, (self.q)] = matrixMult(self.Kve.T, X)[:, 0]
            X = self.X0t_stack

        H = np.array(range(ngrids)) / float(ngrids)
        L = np.array([self.LL(h, X, stack=False, REML=REML)[0] for h in H])
        self.LLs = L

        hmax = self.getMax(H, X, REML)
        L, beta, sigma, betaSTDERR = self.LL(hmax, X, stack=False, REML=REML)

        self.H = H
        self.optH = hmax.sum()
        self.optLL = L
        self.optBeta = beta
        self.optSigma = sigma.sum()

        return hmax, beta, sigma, L

    def association(self, X, h=None, stack=True, REML=True, returnBeta=False, L0=None):

        """
        Calculates association statitics for the SNPs encoded in the vector X of size n.
        If h == None, the optimal h stored in optH is used.
        """
        if stack:
            # X = np.hstack([self.X0t,matrixMult(self.Kve.T, X)])
            self.X0t_stack[:, (self.q)] = matrixMult(self.Kve.T, X)[:, 0]
            X = self.X0t_stack

        if h == None: h = self.optH

        if L0 is not None:
            # likelihood ratio test is only valid for ML estimate
            REML = False

        L, beta, sigma, betaVAR = self.LL(h, X, stack=False, REML=REML)
        q = len(beta)

        if L0 is None:
            ts, ps, ds = self.tstat(beta[q - 1], betaVAR[q - 1, q - 1], sigma, q)

            if returnBeta: return ts, ps, ds, beta[q - 1].sum(), betaVAR[q - 1, q - 1].sum() * sigma
            return ts, ps, ds

        else:
            # perfome log-likelihood ratio test
            LR = -2. * (L0 - L)
            ps = stats.chi2.sf(LR, df=1)
            return LR, ps, np.sqrt(LR / self.N), beta[q - 1].sum(), betaVAR[q - 1, q - 1].sum() * sigma

    def epistasis(self, X, h=None, stack=True):
        """
        Calculates association statistics for a pair of SNPs encoded in the vector X of size n.
        We calculate the full rank model y ~ x1+ x2 + x1*x2
        and test with a likelihood ration test for epistatic association
        """
        if stack:
            # X = np.hstack([self.X0t,matrixMult(self.Kve.T, X)])
            self.X0t_stack[:, (self.q)] = matrixMult(self.Kve.T, X)[:, 0]
            X = self.X0t_stack

        if h is None:
            h = self.optH

        X_0 = X[:, :-1]
        L_0, beta_0, sigma_0, betaVAR_0 = self.LL(h,
                                                  X_0,
                                                  stack=False,
                                                  REML=False)

        L_epi, beta_epi, sigma_epi, betaVAR_epi = self.LL(h,
                                                          X,
                                                          stack=False,
                                                          REML=False)

        # test with a log-likelihood ratio test
        LR = -2. * (L_0 - L_epi)
        ps = stats.chi2.sf(LR, df=1)
        return LR, ps, np.sqrt(LR / self.N), beta_epi[-1]

    def tstat(self, beta, var, sigma, q, log=False):

        """
        Calculates a t-statistic and associated p-value given the estimate of beta and its standard error.
        This is actually an F-test, but when only one hypothesis is being performed, it reduces to a t-test.
        """

        ts = beta / np.sqrt(var * sigma)
        cohens_d = beta * sigma
        # ps = 2.0*(1.0 - stats.t.cdf(np.abs(ts), self.N-q))
        # sf == survival function - this is more accurate -- could also use logsf if the precision is not good enough
        if log:
            ps = 2.0 + (stats.t.logsf(np.abs(ts), self.N - q))
        else:
            ps = 2.0 * (stats.t.sf(np.abs(ts), self.N - q))
        if not len(ts) == 1 or not len(ps) == 1:
            raise Exception("Something bad happened :(")
        return ts.sum(), ps.sum(), cohens_d.sum()

    def plotFit(self, color='b-', title=''):

        """
        Simple function to visualize the likelihood space.  It takes the LLs
        calcualted over a grid and normalizes them by subtracting off the mean and exponentiating.
        The resulting "probabilities" are normalized to one and plotted against heritability.
        This can be seen as an approximation to the posterior distribuiton of heritability.

        For diagnostic purposes this lets you see if there is one distinct maximum or multiple
        and what the variance of the parameter looks like.
        """
        import matplotlib.pyplot as pl

        mx = self.LLs.max()
        p = np.exp(self.LLs - mx)
        p = p / p.sum()

        pl.plot(self.H, p, color)
        pl.xlabel("Heritability")
        pl.ylabel("Probability of data")
        pl.title(title)

    def meanAndVar(self):

        mx = self.LLs.max()
        p = np.exp(self.LLs - mx)
        p = p / p.sum()

        mn = (self.H * p).sum()
        vx = ((self.H - mn) ** 2 * p).sum()

        return mn, vx
