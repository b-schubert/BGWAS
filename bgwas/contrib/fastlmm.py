"""
   Copyright 2014 Microsoft Corporation

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

https://github.com/MicrosoftGenomics/FaST-LMM
"""

import scipy as sp
import numpy as np
import scipy.linalg as la
import scipy.optimize as opt
import scipy.special as ss

import time
import logging


class LMM(object):
    """
    linear mixed model with up to two kernels
    N(y | X*beta ; sigma2(h2*((1-a2)*K0 + a2*K1) + (1-h2)*I),
    where
    K0 = G0*G0^T
    K1 = G1*G1^T
    """
    __slots__ = ["G", "G0", "G1", "y", "X", "K0", "K1", "K", "U", "S", "UX", "Uy", "UUX", "UW", "UUW", "UUy", "pos0",
                 "pos1", "a2", "exclude_idx",
                 "forcefullrank", "numcalls", "Xstar", "Kstar", "Kstar_star", "UKstar", "UUKstar", "Gstar", "K0star",
                 "K1star", "K0star_star", "K1star_star"]

    def __init__(self, forcefullrank=False):
        '''
        Input:
        forcefullrank   : if True, then the code always computes K and runs cubically
                            (False)
        '''
        self.X = None
        self.y = None
        self.G = None
        self.G0 = None
        self.G1 = None
        self.K = None
        self.K0 = None
        self.K1 = None
        self.U = None
        self.S = None
        self.Uy = None
        self.UUy = None
        self.UX = None
        self.UUX = None
        self.UW = None
        self.UUW = None
        self.pos0 = None
        self.pos1 = None
        self.a2 = None
        self.exclude_idx = []
        self.forcefullrank = forcefullrank
        self.numcalls = 0
        self.Xstar = None
        self.Kstar = None
        self.Kstar_star = None
        self.UKstar = None
        self.UUKstar = None
        self.Gstar = None

    def setX(self, X):
        '''
        set the fixed effects X (covariates).
        The Kernel has to be set in advance by first calling setG() or setK().
        --------------------------------------------------------------------------
        Input:
        X       : [N*D] 2-dimensional array of covariates
        --------------------------------------------------------------------------
        '''
        self.X = X
        self.UX = self.U.T.dot(X)
        k = self.S.shape[0]
        N = self.X.shape[0]
        if (k < N):
            self.UUX = X - self.U.dot(self.UX)

    def setX2(self, X):
        '''
        a version of setX that doesn't assume that Eigenvalue decomposition has been done.
        '''
        self.X = X
        N = self.X.shape[0]

    def sety(self, y):
        '''
        set the phenotype y.
        The Kernel has to be set in advance by first calling setG() or setK().
        --------------------------------------------------------------------------
        Input:
        y       : [N] 1-dimensional array of phenotype values
        --------------------------------------------------------------------------
        '''
        assert y.ndim == 1, "y should be 1-dimensional"
        self.y = y
        self.Uy = self.U.T.dot(y)
        k = self.S.shape[0]
        N = self.y.shape[0]
        if (k < N):
            self.UUy = y - self.U.dot(self.Uy)

    def sety2(self, y):
        '''
        a version of sety that doesn't assume that Eigenvalue decomposition has been done.
        '''
        assert y.ndim == 1, "y should be 1-dimensional"
        self.y = y
        N = self.y.shape[0]

    def setG(self, G0=None, G1=None, a2=0.0, K0=None, K1=None):
        '''
        set the Kernel (1-a2)*K0 and a2*K1 from G0 and G1.
        This has to be done before setting the data setX() and setY().

        If k0+k1>>N and similar kernels are used repeatedly, it is beneficial to precompute
        the kernel and pass it as an argument.
        ----------------------------------------------------------------------------
        Input:
        G0              : [N*k0] array of random effects
        G1              : [N*k1] array of random effects (optional)
        a2              : mixture weight between K0=G0*G0^T and K1=G1*G1^T

        K0              : [N*N] array, random effects covariance (positive semi-definite)
        K1              : [N*N] array, random effects covariance (positive semi-definite)(optional)
        -----------------------------------------------------------------------------
        '''
        self.G0 = G0
        self.G1 = G1
        if a2 < 0.0:
            a2 = 0.0
        if a2 > 1.0:
            a2 = 1.0

        if G1 is None and G0 is not None:
            self.G = G0
        elif G0 is not None and G1 is not None:
            # build the weighted concatenation of G0 and G1 = varianceComponent
            if a2 == 0.0:
                logging.info("a2=0.0, only using G0")
                self.G = G0
            elif a2 == 1.0:
                self.G = G1
                logging.info("a2=1.0, only using G1")
            else:
                self.G = sp.concatenate((sp.sqrt(1.0 - a2) * G0, sp.sqrt(a2) * G1), 1)

        else:
            self.G = None

        if self.G is not None:
            N = self.G.shape[0]
            k = self.G.shape[1]
        else:
            N = K0.shape[0]
            k = N
        if k > 0:
            if ((not self.forcefullrank) and (k < N)):
                # it is faster using the eigen decomposition of G.T*G but this is more accurate
                try:
                    [U, S, V] = la.svd(self.G, full_matrices=False)
                    if np.any(S < -0.1):
                        logging.warning("kernel contains a negative Eigenvalue")
                    self.U = U
                    self.S = S * S

                except la.LinAlgError:  # revert to Eigenvalue decomposition
                    logging.warning(
                        "Got SVD exception, trying eigenvalue decomposition of square of G. \
                         Note that this is a little bit less accurate")
                    [S_, V_] = la.eigh(self.G.T.dot(self.G))
                    if np.any(S_ < -0.1):
                        logging.warning("kernel contains a negative Eigenvalue")
                    S_nonz = (S_ > 0)
                    self.S = S_[S_nonz]
                    self.S *= (N / self.S.sum())
                    self.U = self.G.dot(V_[:, S_nonz] / sp.sqrt(self.S))
            else:
                if K0 is None:
                    K0 = self.G0.dot(self.G0.T)
                self.K0 = K0
                if (self.G1 is not None) and (K1 is None):
                    K1 = self.G1.dot(self.G1.T)
                self.setK(K0=K0, K1=K1, a2=a2)
                # K=self.G.dot(self.G.T)
                # self.setK(K)
            self.a2 = a2
            pass
        else:  # rank of kernel = 0 (linear regression case)
            self.S = sp.zeros((0))
            self.U = sp.zeros_like(self.G)

    def setK(self, K0, K1=None, a2=0.0):
        '''
        set the Kernel (1-a2)*K0 and a2*K1.
        This has to be done before setting the data setX() and setY().
        --------------------------------------------------------------------------
        Input:
        K0 : [N*N] array, random effects covariance (positive semi-definite)
        K1 : [N*N] array, random effects covariance (positive semi-definite)(optional)
        a2 : mixture weight between K0 and K1
        --------------------------------------------------------------------------
        '''
        self.K0 = K0
        self.K1 = K1
        logging.debug("About to mix K0 and K1")
        if K1 is None:
            self.K = K0
        else:
            self.K = (1.0 - a2) * K0 + a2 * K1
        logging.debug("About to eigh")
        [S, U] = la.eigh(self.K)
        logging.debug("Done with to eigh")
        if np.any(S < -0.1):
            logging.warning("kernel contains a negative Eigenvalue")

        self.U = U
        self.S = S  # *(S.shape[0]/S.sum())
        self.a2 = a2

    def setK2(self, K0, K1=None, a2=0.0):
        '''
        a version of setK that doesn't do Eigenvalue decomposition.
        '''
        self.K0 = K0
        self.K1 = K1
        logging.debug("About to mix K0 and K1")
        if K1 is None:
            self.K = K0
        else:
            self.K = (1.0 - a2) * K0 + a2 * K1
        self.a2 = a2

    def set_exclude_idx(self, idx):
        '''

        --------------------------------------------------------------------------
        Input:
        idx  : [k_up: number of SNPs to be removed] holds the indices of SNPs to be removed
        --------------------------------------------------------------------------
        '''

        self.exclude_idx = idx

    def innerLoopTwoKernel(self, a2=0.5, nGridH2=10, minH2=0.0, maxH2=0.99999, **kwargs):
        '''
        For a given weight a2, finds the optimal h2 and returns the negative log-likelihood
        --------------------------------------------------------------------------
        Input:
        a2      : mixture weight between K0 and K1
        nGridH2 : number of h2-grid points to evaluate the negative log-likelihood at
        minH2   : minimum value for h2 optimization
        maxH2   : maximum value for h2 optimization
        --------------------------------------------------------------------------
        Output:
        dictionary containing the model parameters at the optimal h2
        --------------------------------------------------------------------------
        '''

        if self.K0 is not None:
            self.setK(K0=self.K0, K1=self.K1, a2=a2)
        else:
            self.setG(G0=self.G0, G1=self.G1, a2=a2)
        self.setX(self.X)
        self.sety(self.y)
        return self.findH2(nGridH2=nGridH2, minH2=minH2, maxH2=maxH2, **kwargs)

    def findA2(self, nGridA2=10, minA2=0.0, maxA2=1.0, nGridH2=10, minH2=0.0, maxH2=0.99999, verbose=False, **kwargs):
        '''
        Find the optimal a2 and h2, such that K=(1.0-a2)*K0+a2*K1. Performs a double loop optimization (could be expensive for large grid-sizes)
        (default maxA2 value is set to 1 as loss of positive definiteness of the final model covariance only depends on h2, not a2)
        --------------------------------------------------------------------------
        Input:
        nGridA2 : number of a2-grid points to evaluate the negative log-likelihood at
        minA2   : minimum value for a2 optimization
        maxA2   : maximum value for a2 optimization
        nGridH2 : number of h2-grid points to evaluate the negative log-likelihood at
        minH2   : minimum value for h2 optimization
        maxH2   : maximum value for h2 optimization
        --------------------------------------------------------------------------
        Output:
        dictionary containing the model parameters at the optimal h2 and a2
        --------------------------------------------------------------------------
        '''
        self.numcalls = 0
        resmin = [None]

        def f(x, resmin=resmin, nGridH2=nGridH2, minH2=minH2, maxH2=maxH2, **kwargs):
            self.numcalls += 1
            t0 = time.time()
            res = self.innerLoopTwoKernel(a2=x, nGridH2=nGridH2, minH2=minH2, maxH2=maxH2, **kwargs)
            if (resmin[0] is None) or (res['nLL'] < resmin[0]['nLL']):
                resmin[0] = res
            t1 = time.time()
            logging.info("x={0}. one objective function call took {1} seconds elapsed ".format(x, t1 - t0))
            # import pdb; pdb.set_trace()
            return res['nLL']

        if verbose: logging.info("finda2")
        min = minimize1D(f=f, nGrid=nGridA2, minval=minA2, maxval=maxA2, verbose=False)
        # print "numcalls to innerLoopTwoKernel= " + str(self.numcalls)
        return resmin[0]

    def findH2(self, nGridH2=10, minH2=0.0, maxH2=0.99999, **kwargs):
        '''
        Find the optimal h2 for a given K. Note that this is the single kernel case. So there is no a2.
        (default maxH2 value is set to a value smaller than 1 to avoid loss of positive definiteness of the final model covariance)
        --------------------------------------------------------------------------
        Input:
        nGridH2 : number of h2-grid points to evaluate the negative log-likelihood at
        minH2   : minimum value for h2 optimization
        maxH2   : maximum value for h2 optimization
        --------------------------------------------------------------------------
        Output:
        dictionary containing the model parameters at the optimal h2
        --------------------------------------------------------------------------
        '''
        # f = lambda x : (self.nLLeval(h2=x,**kwargs)['nLL'])
        resmin = [None]

        def f(x, resmin=resmin, **kwargs):
            res = self.nLLeval(h2=x, **kwargs)
            if (resmin[0] is None) or (res['nLL'] < resmin[0]['nLL']):
                resmin[0] = res
            logging.debug("search\t{0}\t{1}".format(x, res['nLL']))
            return res['nLL']

        min = minimize1D(f=f, nGrid=nGridH2, minval=minH2, maxval=maxH2)
        return resmin[0]

    def find_log_delta(self, sid_count, min_log_delta=-5, max_log_delta=10, nGrid=100, **kwargs):
        '''
        #Need comments
        '''
        # f = lambda x : (self.nLLeval(h2=x,**kwargs)['nLL'])
        resmin = [None]

        def f(x, resmin=resmin, **kwargs):
            # We convert from external log_delta to h2 and then back again so that this code is most similar to findH2
            h2 = 1.0 / (np.exp(x) * sid_count + 1)
            res = self.nLLeval(h2=h2, **kwargs)
            if (resmin[0] is None) or (res['nLL'] < resmin[0]['nLL']):
                resmin[0] = res
            return res['nLL']

        min = minimize1D(f=f, nGrid=nGrid, minval=min_log_delta, maxval=max_log_delta)
        res = resmin[0]
        internal_delta = 1.0 / res['h2'] - 1.0
        ln_external_delta = np.log(internal_delta / sid_count)
        res['log_delta'] = ln_external_delta
        return res

    def nLLeval(self, h2=0.0, REML=True, logdelta=None, delta=None, dof=None, scale=1.0, penalty=0.0):
        '''
        evaluate -ln( N( U^T*y | U^T*X*beta , h2*S + (1-h2)*I ) ),
        where ((1-a2)*K0 + a2*K1) = USU^T
        --------------------------------------------------------------------------
        Input:
        h2      : mixture weight between K and Identity (environmental noise)
        REML    : boolean
                  if True   : compute REML
                  if False  : compute ML
        dof     : Degrees of freedom of the Multivariate student-t
                        (default None uses multivariate Normal likelihood)
        logdelta: log(delta) allows to optionally parameterize in delta space
        delta   : delta     allows to optionally parameterize in delta space
        scale   : Scale parameter the multiplies the Covariance matrix (default 1.0)
        --------------------------------------------------------------------------
        Output dictionary:
        'nLL'       : negative log-likelihood
        'sigma2'    : the model variance sigma^2
        'beta'      : [D*1] array of fixed effects weights beta
        'h2'        : mixture weight between Covariance and noise
        'REML'      : True: REML was computed, False: ML was computed
        'a2'        : mixture weight between K0 and K1
        'dof'       : Degrees of freedom of the Multivariate student-t
                        (default None uses multivariate Normal likelihood)
        'scale'     : Scale parameter that multiplies the Covariance matrix (default 1.0)
        --------------------------------------------------------------------------
        '''
        if (h2 < 0.0) or (h2 > 1.0):
            return {'nLL': 3E20,
                    'h2': h2,
                    'REML': REML,
                    'scale': scale}
        k = self.S.shape[0]
        N = self.y.shape[0]
        D = self.UX.shape[1]

        # if REML == True:
        #    # this needs to be fixed, please see test_gwas.py for details
        #    raise NotImplementedError("this feature is not ready to use at this time, please use lmm_cov.py instead")

        if logdelta is not None:
            delta = sp.exp(logdelta)

        if delta is not None:
            Sd = (self.S + delta) * scale
        else:
            Sd = (h2 * self.S + (1.0 - h2)) * scale

        UXS = self.UX / np.lib.stride_tricks.as_strided(Sd, (Sd.size, self.UX.shape[1]), (Sd.itemsize, 0))
        UyS = self.Uy / Sd

        XKX = UXS.T.dot(self.UX)
        XKy = UXS.T.dot(self.Uy)
        yKy = UyS.T.dot(self.Uy)

        logdetK = sp.log(Sd).sum()

        if (k < N):  # low rank part

            # determine normalization factor
            if delta is not None:
                denom = (delta * scale)
            else:
                denom = ((1.0 - h2) * scale)

            XKX += self.UUX.T.dot(self.UUX) / (denom)
            XKy += self.UUX.T.dot(self.UUy) / (denom)
            yKy += self.UUy.T.dot(self.UUy) / (denom)
            logdetK += (N - k) * sp.log(denom)

        # proximal contamination (see Supplement Note 2: An Efficient Algorithm for Avoiding Proximal Contamination)
        # available at: http://www.nature.com/nmeth/journal/v9/n6/extref/nmeth.2037-S1.pdf
        # exclude SNPs from the RRM in the likelihood evaluation

        if len(self.exclude_idx) > 0:
            num_exclude = len(self.exclude_idx)

            # consider only excluded SNPs
            G_exclude = self.G[:, self.exclude_idx]

            self.UW = self.U.T.dot(G_exclude)  # needed for proximal contamination
            UWS = self.UW / np.lib.stride_tricks.as_strided(Sd, (Sd.size, num_exclude), (Sd.itemsize, 0))
            assert UWS.shape == (k, num_exclude)

            WW = np.eye(num_exclude) - UWS.T.dot(self.UW)
            WX = UWS.T.dot(self.UX)
            Wy = UWS.T.dot(self.Uy)
            assert WW.shape == (num_exclude, num_exclude)
            assert WX.shape == (num_exclude, D)
            assert Wy.shape == (num_exclude,)

            if (k < N):  # low rank part

                self.UUW = G_exclude - self.U.dot(self.UW)

                WW += self.UUW.T.dot(self.UUW) / denom
                WX += self.UUW.T.dot(self.UUX) / denom
                Wy += self.UUW.T.dot(self.UUy) / denom

            # TODO: do cholesky, if fails do eigh
            # compute inverse efficiently
            [S_WW, U_WW] = la.eigh(WW)

            UWX = U_WW.T.dot(WX)
            UWy = U_WW.T.dot(Wy)
            assert UWX.shape == (num_exclude, D)
            assert UWy.shape == (num_exclude,)

            # compute S_WW^{-1} * UWX
            WX = UWX / np.lib.stride_tricks.as_strided(S_WW, (S_WW.size, UWX.shape[1]), (S_WW.itemsize, 0))
            # compute S_WW^{-1} * UWy
            Wy = UWy / S_WW
            # determinant update
            logdetK += sp.log(S_WW).sum()
            assert WX.shape == (num_exclude, D)
            assert Wy.shape == (num_exclude,)

            # perform updates (instantiations for a and b in Equation (1.5) of Supplement)
            yKy += UWy.T.dot(Wy)
            XKy += UWX.T.dot(Wy)
            XKX += UWX.T.dot(WX)

        #######

        [SxKx, UxKx] = la.eigh(XKX)
        # optionally regularize the beta weights by penalty
        if penalty > 0.0:
            SxKx += penalty
        i_pos = SxKx > 1E-10
        beta = sp.dot(UxKx[:, i_pos], (sp.dot(UxKx[:, i_pos].T, XKy) / SxKx[i_pos]))

        r2 = yKy - XKy.dot(beta)

        if dof is None:  # Use the Multivariate Gaussian
            if REML:
                XX = self.X.T.dot(self.X)
                [Sxx, Uxx] = la.eigh(XX)
                logdetXX = sp.log(Sxx).sum()
                logdetXKX = sp.log(SxKx).sum()
                sigma2 = r2 / (N - D)
                nLL = 0.5 * (logdetK + logdetXKX - logdetXX + (N - D) * (sp.log(2.0 * sp.pi * sigma2) + 1))
            else:
                sigma2 = r2 / (N)
                nLL = 0.5 * (logdetK + N * (sp.log(2.0 * sp.pi * sigma2) + 1))
            result = {
                'nLL': nLL,
                'sigma2': sigma2,
                'beta': beta,
                'h2': h2,
                'REML': REML,
                'a2': self.a2,
                'scale': scale
            }
        else:  # Use multivariate student-t
            if REML:
                XX = self.X.T.dot(self.X)
                [Sxx, Uxx] = la.eigh(XX)
                logdetXX = sp.log(Sxx).sum()
                logdetXKX = sp.log(SxKx).sum()

                nLL = 0.5 * (logdetK + logdetXKX - logdetXX + (dof + (N - D)) * sp.log(1.0 + r2 / dof))
                nLL += 0.5 * (N - D) * sp.log(dof * sp.pi) + ss.gammaln(0.5 * dof) - ss.gammaln(0.5 * (dof + (N - D)))
            else:
                nLL = 0.5 * (logdetK + (dof + N) * sp.log(1.0 + r2 / dof))
                nLL += 0.5 * N * sp.log(dof * sp.pi) + ss.gammaln(0.5 * dof) - ss.gammaln(0.5 * (dof + N))
            result = {
                'nLL': nLL,
                'dof': dof,
                'beta': beta,
                'h2': h2,
                'REML': REML,
                'a2': self.a2,
                'scale': scale
            }
        assert sp.all(sp.isreal(nLL)), "nLL has an imaginary component, possibly due to constant covariates"
        return result

    def getPosteriorWeights(self, beta, h2=0.0, logdelta=None, delta=None, scale=1.0):
        '''
        compute posterior mean over the feature weights (effect sizes of SNPs in the kernel, not the SNPs being tested):
        w = G.T (GG.T + delta*I)^(-1) (y - Xbeta)
        --------------------------------------------------------------------------
        Input:
        beta            : weight vector for fixed effects
        h2              : mixture weight between K and Identity (environmental noise)
        logdelta        : log(delta) allows to optionally parameterize in delta space
        delta           : delta     allows to optionally parameterize in delta space
        scale           : Scale parameter the multiplies the Covariance matrix (default 1.0)

        returnVar       : if True, marginal variances are estimated
        returnCovar     : if True, posterior covariance is learnt
        --------------------------------------------------------------------------
        Dictionary with the following fields:
        weights           : [k0+k1] 1-dimensional array of predicted phenotype values
        --------------------------------------------------------------------------
        '''
        k = self.S.shape[0]
        N = self.y.shape[0]

        if logdelta is not None:
            delta = sp.exp(logdelta)
        if delta is not None:
            Sd = (self.S + delta) * scale
        else:
            Sd = (h2 * self.S + (1.0 - h2)) * scale

        yres = self.y - sp.dot(self.X, beta)
        Uyres = sp.dot(self.U.T, yres)
        UG = sp.dot(self.U.T, self.G)
        weights = sp.dot(UG.T, Uyres / Sd)

        if k < N:  # low-rank part
            # determine normalization factor
            if delta is not None:
                denom = (delta * scale)
            else:
                denom = ((1.0 - h2) * scale)

            UUG = self.G - sp.dot(self.U, UG)
            UUyres = yres - sp.dot(self.U, Uyres)
            weights += UUG.T.dot(UUyres) / (denom)

        return weights

    def setTestData(self, Xstar, K0star=None, K1star=None, G0star=None, G1star=None):
        '''
        set data for predicting

        --------------------------------------------------------------------------
        Input:
        Xstar           : [M,D] 2-dimensional array of covariates on the test set
        G0star          : [M,k0] array of random effects on the test set
        G1star          : [M,k1] array of random effects on the test set (optional)
        K0star          : [M,N] array, random effects covariance between test and training data (positive semi-definite)
        K1star          : [M,N] array, random effects covariance between test and training data (positive semi-definite)(optional)
        where M is # of test cases, N is the # of training cases
        --------------------------------------------------------------------------
        '''

        self.Xstar = Xstar
        if G1star is None:
            self.Gstar = G0star
        else:
            if self.a2 == 0.0:
                logging.info("a2=0.0, only using G0")
                self.Gstar = G0star
            elif self.a2 == 1.0:
                self.Gstar = G1star
                logging.info("a2=1.0, only using G1")
            else:
                self.Gstar = sp.concatenate((sp.sqrt(1.0 - self.a2) * G0star, sp.sqrt(self.a2) * G1star), 1)

        if K0star is not None:
            if K1star is None:
                self.Kstar = K0star
            else:
                self.Kstar = (1.0 - self.a2) * K0star + self.a2 * K1star
        else:
            self.Kstar = sp.dot(self.Gstar, self.G.T)

        self.UKstar = sp.dot(self.U.T, self.Kstar.T)

        if self.G is not None:
            k = self.G.shape[1]
            N = self.G.shape[0]
            if k < N:
                # see e.g. Equation 3.17 in Supplement of FaST LMM paper
                self.UUKstar = self.Kstar.T - sp.dot(self.U, self.UKstar)

    def setTestData2(self, Xstar, K0star=None, K1star=None):
        '''
        a version of setTestData that doesn't assume that Eigenvalue decomposition has been done.
        '''

        self.Xstar = Xstar
        self.Gstar = None
        if K1star is None:
            self.Kstar = K0star
        else:
            self.Kstar = (1.0 - self.a2) * K0star + self.a2 * K1star

    def predictMean(self, beta, h2=0.0, logdelta=None, delta=None, scale=1.0):
        '''
        mean prediction for the linear mixed model on unobserved data:

        ystar = X*beta + Kstar(h2*K + (1-h2)*K)^{-1}(y-X*beta)
        where Kstar is the train vs test kernel
        --------------------------------------------------------------------------
        Input:
        beta            : weight vector for fixed effects
        h2              : mixture weight between K and Identity (environmental noise)
        logdelta        : log(delta) allows to optionally parameterize in delta space
        delta           : delta     allows to optionally parameterize in delta space
        scale           : Scale parameter the multiplies the Covariance matrix (default 1.0)


        If SNPs are excluded, nLLeval must be called before to re-calculate self.UW,self.UUW
        --------------------------------------------------------------------------
        Output:
        ystar           : [M] 1-dimensional array of predicted phenotype values
        --------------------------------------------------------------------------
        '''

        M = self.Xstar.shape[0]

        if (h2 < 0.0) or (h2 >= 1.0):
            return sp.nan * sp.ones(M)

        k = self.S.shape[0]
        N = self.y.shape[0]
        # D=self.UX.shape[1]

        if logdelta is not None:
            delta = sp.exp(logdelta)

        # delta = (1-h2) / h2
        if delta is not None:
            Sd = (self.S + delta) * scale
        else:
            assert False, "not implemented (UKstar needs to be scaled by h2)"
            Sd = (h2 * self.S + (1.0 - h2)) * scale

        if len(self.exclude_idx) > 0:
            # cut out
            num_exclude = len(self.exclude_idx)
            # consider only excluded SNPs
            Gstar_exclude = self.Gstar[:, self.exclude_idx]
            # G_exclude = self.G[:,self.exclude_idx]
            UKstar = self.UKstar - sp.dot(self.UW, Gstar_exclude.T)
            if k < N:
                UUKstar = self.UUKstar - sp.dot(self.UUW, Gstar_exclude.T)
        else:
            UKstar = self.UKstar
            UUKstar = self.UUKstar

        yfixed = sp.dot(self.Xstar, beta)
        yres = self.y - sp.dot(self.X, beta)
        Uyres = self.Uy - sp.dot(self.UX, beta)
        Sdi = 1. / Sd
        yrandom = sp.dot(Sdi * UKstar.T, Uyres)

        if k < N:  # low-rank part
            # determine normalization factor
            if delta is not None:
                denom = (delta * scale)
            else:
                denom = ((1.0 - h2) * scale)
            UUyres = yres - sp.dot(self.U, Uyres)
            yrandom += sp.dot(UUKstar.T, UUyres) / denom

        # proximal contamination (see Supplement Note 2: An Efficient Algorithm for Avoiding Proximal Contamination)
        # available at: http://www.nature.com/nmeth/journal/v9/n6/extref/nmeth.2037-S1.pdf
        # exclude SNPs from the RRM in the likelihood evaluation
        if len(self.exclude_idx) > 0:
            UWS = self.UW / np.lib.stride_tricks.as_strided(Sd, (Sd.size, num_exclude), (Sd.itemsize, 0))
            assert UWS.shape == (k, num_exclude)
            WW = np.eye(num_exclude) - UWS.T.dot(self.UW)
            WKstar = UWS.T.dot(UKstar)
            Wyres = UWS.T.dot(Uyres)
            assert WW.shape == (num_exclude, num_exclude)
            assert WKstar.shape == (num_exclude, M)
            assert Wyres.shape == (num_exclude,)

            if (k < N):  # low rank part
                WW += self.UUW.T.dot(self.UUW) / denom
                WKstar += self.UUW.T.dot(UUKstar) / denom
                Wyres += self.UUW.T.dot(UUyres) / denom

            # TODO: do cholesky, if fails do eigh
            # compute inverse efficiently
            [S_WW, U_WW] = la.eigh(WW)

            UWKstar = U_WW.T.dot(WKstar)
            UWyres = U_WW.T.dot(Wyres)
            assert UWKstar.shape == (num_exclude, M)
            assert UWyres.shape == (num_exclude,)

            # compute S_WW^{-1} * UWX
            WKstar = UWKstar / np.lib.stride_tricks.as_strided(S_WW, (S_WW.size, UWKstar.shape[1]), (S_WW.itemsize, 0))
            # compute S_WW^{-1} * UWy
            Wyres = UWyres / S_WW
            assert WKstar.shape == (num_exclude, M)
            assert Wyres.shape == (num_exclude,)

            # perform updates (instantiations for a and b in Equation (1.5) of Supplement)
            yrandom += UWKstar.T.dot(Wyres)

        ystar = yfixed + yrandom
        return ystar

    def predict_mean_and_variance(lmm, beta, sigma2, h2, Kstar_star):
        assert 0 <= h2 and h2 <= 1, "By definition, h2 must be between 0 and 1 (inclusive)"
        varg = h2 * sigma2
        vare = (1. - h2) * sigma2
        if lmm.G is not None:
            K = np.dot(lmm.G, lmm.G.T)  # !!!later this is very inefficient in memory and computation
        else:
            K = np.dot(np.dot(lmm.U, np.eye(len(lmm.U)) * lmm.S),
                       lmm.U.T)  # Re-compose the Eigen value decomposition #!!!later do this more efficiently
        V = varg * K + vare * np.eye(len(K))
        Vinv = la.inv(V)

        a = np.dot(varg * lmm.Kstar, Vinv)

        y_star = np.dot(lmm.Xstar, beta) + np.dot(a, lmm.y - sp.dot(lmm.X,
                                                                    beta))  # !!!later shouldn't the 2nd dot be precomputed?
        y_star = y_star.reshape(-1, 1)  # Make 2-d

        var_star = (varg * Kstar_star +
                    vare * np.eye(len(Kstar_star)) -
                    np.dot(a,
                           (varg * lmm.Kstar.T)))
        return y_star, var_star

    def nLL(lmm, beta, sigma2, h2, y_actual):
        from scipy.stats import multivariate_normal
        # TODO: predict_mean_and_variance
        y_star, var_star = predict_mean_and_variance(lmm, beta, sigma2, h2, lmm.Kstar_star)
        var = multivariate_normal(mean=y_star.reshape(-1), cov=var_star)
        return -np.log(var.pdf(y_actual.reshape(-1)))

    def predictVariance(self, h2=0.0, logdelta=None, delta=None, sigma2=1.0, Kstar_star=None):
        '''
        variance prediction for the linear mixed model on unobserved data:

        Var_star = sigma2 * (K(X*,X*) + delta*I - Kstar (K + delta*I)^{-1} Kstar )
        --------------------------------------------------------------------------
        Input:
        h2              : mixture weight between K and Identity (environmental noise)
        logdelta        : log(delta) allows to optionally parameterize in delta space
        delta           : delta     allows to optionally parameterize in delta space
        sigma2          : sigma2 parameter the multiplies the Covariance matrix (default 1.0)
        K_star_star     : Kernel on test examples

        If SNPs are excluded, nLLeval must be called before to re-calculate self.UW,self.UUW
        --------------------------------------------------------------------------
        Output:
        Cov_star           : [M,M] 2-dimensional array covariance matrix
        --------------------------------------------------------------------------
        '''

        # TODO: proximal contamination
        # TODO: REML?

        if (h2 < 0.0) or (h2 >= 1.0):
            return sp.nan * sp.ones(M)

        k = self.S.shape[0]
        N = self.y.shape[0]
        # D = self.UX.shape[1]

        # print "k, N, D", k, N, D

        if logdelta is not None:
            delta = sp.exp(logdelta)

        if delta is not None:
            # Sd = (self.S+delta)*sigma2
            Sd = (self.S + delta)
        else:
            # Sd = (h2*self.S + (1.0-h2))*sigma2
            Sd = (h2 * self.S + (1.0 - h2))
            assert False, "h2 code path not test. Please use delta or logdelta"
            # delta = 1.0/h2-1.0 #right?

        Sdi = 1. / Sd

        # part 1 from c-code
        # TODO: handle h2 parameterization
        # TODO: make more efficient (add_diag)

        if Kstar_star is None:
            N_test = self.Gstar.shape[0]
            Kstar_star = sp.dot(self.Gstar, self.Gstar.T)
        else:
            Kstar_star = Kstar_star.copy()

        N_test = Kstar_star.shape[0]
        assert N_test == Kstar_star.shape[1]

        part1 = Kstar_star
        part1 += sp.eye(N_test) * delta
        part1 *= sigma2

        # print "part1", part1[0,0]
        # print "delta", delta, "sigma2", sigma2

        # part 2 from c-code
        # (U1^T a)^T (S_1 + delta*I)^{-1} (U1^T a)
        SUKstarTUkStar = sp.dot(Sdi * self.UKstar.T, self.UKstar)

        # UXS = self.UKstar / NP.lib.stride_tricks.as_strided(Sd, (Sd.size,self.UKstar.shape[1]), (Sd.itemsize,0))
        # NP.testing.assert_array_almost_equal(SUKstarTUkStar, SP.dot(UXS.T, self.UKstar), decimal=4)

        SUKstarTUkStar *= sigma2

        # print "UKstar[0,0]", self.UKstar[0,0]
        # print "UKstarS[0,0]", UXS[0,0]
        # print "SUK", SUKstarTUkStar[0,0]

        # part 3&4 from c-code
        if k < N:  # low-rank part
            # determine normalization factor
            if delta is not None:
                denom = (delta * sigma2)
            else:
                denom = ((1.0 - h2) * sigma2)

            # see Equation 3.17 in Supplement of FaST LMM paper:
            # 1 / delta * (((I_n - U1U1^T)a)^T (I_n - U1U1^T)a), a=K(XS,X)
            SUKstarTUkStar += sp.dot(self.UUKstar.T, self.UUKstar) / denom

        # see Carl Rasmussen's book on GPs, Equation 2.24
        # or Equation 5 in Lasso-LMM paper
        Var_star = part1 - SUKstarTUkStar

        return Var_star

    def nLLeval_test(self, y_test, beta, h2=0.0, logdelta=None, delta=None, sigma2=1.0, Kstar_star=None, robust=False):
        """
        compute out-of-sample log-likelihood

        robust: boolean
                indicates if eigenvalues will be truncated at 1E-9 or 1E-4. The former (default) one was used in FastLMMC,
                but may lead to numerically unstable solutions.
        """
        assert y_test.ndim == 1, "y_test should have 1 dimension"
        mu = self.predictMean(beta, h2=h2, logdelta=logdelta, delta=delta)
        res = y_test - mu

        sigma = self.predictVariance(h2=h2, logdelta=logdelta, delta=delta, sigma2=sigma2, Kstar_star=Kstar_star)

        # TODO: benchmark, record speed difference
        """
        # efficient computation of: (y - mu)^T sigma2^{-1} (y - mu)
        # Solve the linear system x = (L L^T)^-1 res

        try:
            L = SP.linalg.cho_factor(sigma)
            res_sig = SP.linalg.cho_solve(L, res)
            logdetK = NP.linalg.slogdet(sigma)[1]

        except Exception, detail:
            print "Cholesky failed, using eigen-value decomposition!"
        """

        [S_, U_] = la.eigh(sigma)

        if robust:
            S_nonz = (S_ > 1E-4)
        else:
            S_nonz = (S_ > 1E-9)
        assert sum(S_nonz) > 0, "Some eigenvalues should be nonzero"
        S = S_[S_nonz]
        U = U_[:, S_nonz]
        Sdi = 1 / S

        res_sig = res.T.dot(Sdi * U).dot(U.T)
        logdetK = sp.log(S).sum()

        # some sanity checks
        if False:
            res_sig3 = sp.linalg.pinv(sigma).dot(res)
            np.testing.assert_array_almost_equal(res_sig, res_sig3, decimal=2)

        # see Carl Rasmussen's book on GPs, equation 5.10, or
        term1 = -0.5 * logdetK
        term2 = -0.5 * sp.dot(res_sig.reshape(-1).T, res.reshape(
            -1))  # Change the inputs to the functions so that these are vectors, not 1xn,nx1
        term3 = -0.5 * len(res) * sp.log(2 * sp.pi)

        if term2 < -10000:
            logging.warning("looks like nLLeval_test is running into numerical difficulties")

            SC = S.copy()
            SC.sort()

            logging.warning(["delta:", delta, "log det", logdetK, "term 2", term2, "term 3:", term3])
            logging.warning(["largest eigv:", SC[-1], "second largest eigv:", SC[-2], "smallest eigv:", SC[0]])
            logging.warning(["ratio 1large/2large:", SC[-1] / SC[-2], "ratio lrg/small:", SC[-1] / SC[0]])

        neg_log_likelihood = -(term1 + term2 + term3)

        return neg_log_likelihood


def minimize1D(f, evalgrid=None, nGrid=10, minval=0.0, maxval=0.99999, verbose=False, brent=True, check_boundaries=True,
               resultgrid=None, return_grid=False):
    '''
    minimize a function f(x) in the grid between minval and maxval.
    The function will be evaluated on a grid and then all triplets,
    where the inner value is smaller than the two outer values are optimized by
    Brent's algorithm.
    --------------------------------------------------------------------------
    Input:
    f(x)    : callable target function
    evalgrid: 1-D array prespecified grid of x-values
    nGrid   : number of x-grid points to evaluate f(x)
    minval  : minimum x-value for optimization of f(x)
    maxval  : maximum x-value for optimization of f(x)
    brent   : boolean indicator whether to do Brent search or not.
              (default: True)
    --------------------------------------------------------------------------
    Output list:
    [xopt, f(xopt)]
    xopt    : x-value at the optimum
    f(xopt) : function value at the optimum
    --------------------------------------------------------------------------
    '''

    if evalgrid is not None and brent:  # if brent we need to sort the input values
        i_sort = evalgrid.argsort()
        evalgrid = evalgrid[i_sort]
    if resultgrid is None:
        [evalgrid, resultgrid] = evalgrid1D(f, evalgrid=evalgrid, nGrid=nGrid, minval=minval, maxval=maxval)

    i_currentmin = resultgrid.argmin()
    minglobal = (evalgrid[i_currentmin], resultgrid[i_currentmin])
    if brent:  # do Brent search in addition to rest?

        if check_boundaries:
            # if the outer boundary point is a local optimum expand search bounded between the grid points
            if resultgrid[0] < resultgrid[1]:
                minlocal = opt.fminbound(f, evalgrid[0], evalgrid[1], full_output=True)
                if minlocal[1] < minglobal[1]:
                    minglobal = minlocal[0:2]

            # if the outer boundary point is a local optimum expand search bounded between the grid points
            if resultgrid[-1] < resultgrid[-2]:
                minlocal = opt.fminbound(f, evalgrid[-2], evalgrid[-1], full_output=True)
                if minlocal[1] < minglobal[1]:
                    minglobal = minlocal[0:2]

        onebrent = False
        # if any triplet is found, where the inner point is a local optimum expand search
        for i in range(resultgrid.shape[0] - 2):
            if (resultgrid[i + 1] < resultgrid[i + 2]) and (resultgrid[i + 1] < resultgrid[i]):
                onebrent = True
                minlocal = opt.brent(f, brack=(evalgrid[i], evalgrid[i + 1], evalgrid[i + 2]), full_output=True)

                if minlocal[1] < minglobal[1]:
                    minglobal = minlocal[0:2]

    if return_grid:
        return minglobal[0], minglobal[1], evalgrid, resultgrid
    else:
        return minglobal


def evalgrid1D(f, evalgrid=None, nGrid=10, minval=0.0, maxval=0.99999, dimF=0):
    '''
    evaluate a function f(x) on all values of a grid.
    --------------------------------------------------------------------------
    Input:
    f(x)    : callable target function
    evalgrid: 1-D array prespecified grid of x-values
    nGrid   : number of x-grid points to evaluate f(x)
    minval  : minimum x-value for optimization of f(x)
    maxval  : maximum x-value for optimization of f(x)
    --------------------------------------------------------------------------
    Output:
    evalgrid    : x-values
    resultgrid  : f(x)-values
    --------------------------------------------------------------------------
    '''
    if evalgrid is None:
        step = (maxval - minval) / (nGrid)
        evalgrid = sp.arange(minval, maxval + step, step)
    if dimF:
        resultgrid = sp.ones((evalgrid.shape[0], dimF)) * 9999999999999.0
    else:
        resultgrid = sp.ones(evalgrid.shape[0]) * 9999999999999.0
    for i in range(evalgrid.shape[0]):
        fevalgrid = f(evalgrid[i])
        assert sp.isreal(fevalgrid).all(), "function returned imaginary value"
        resultgrid[i] = fevalgrid
    return evalgrid, resultgrid
