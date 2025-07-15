import multiprocessing
import numpy as np
import scipy.interpolate
import scipy.integrate
import scipy.special
import emcee
import os

os.environ["OMP_NUM_THREADS"] = "1"

def schechter_function(x=None, alpha=None):
    """Schechter function with unit scaling"""
    f = x**(- alpha) * np.exp(- x)
    return f


def Fabove(alpha=None, log_lambda_star=None, log_lambda_min=None, log_lambda_max=1.,
           log_lambda_c=-2.5):
    integral = scipy.integrate.quad(
        lambda y: 10.**y * np.log(10.) * schechter_function(x=10.**(y - log_lambda_star), alpha=alpha),
        log_lambda_min, log_lambda_max)[0]
    integral_above = scipy.integrate.quad(
        lambda y: 10.**y * np.log(10.) * schechter_function(x=10.**(y - log_lambda_star), alpha=alpha),
        log_lambda_c, log_lambda_max)[0]
    return(integral_above / integral)


def power_law_sample(size=None, alpha=None, log_x_min=None, log_x_max=None):
    r = np.random.random(size)
    map1 = - alpha + 1.
    x_min_term = 10.**(map1 * log_x_min)
    x_max_term = 10.**(map1 * log_x_max)
    x = ((x_max_term - x_min_term) * r + x_min_term)**(1. / map1)
    return(x)


def schechter_sample(size=None, alpha=None, log_x_min=None, log_x_max=None):
    nkeep = 0
    xschechter = np.zeros(size, dtype=np.float32)
    while(nkeep < size):
        ntry = size - nkeep
        xpl = power_law_sample(size=ntry, alpha=alpha, log_x_min=log_x_min,
                               log_x_max=log_x_max)
        choose = np.random.random(ntry)
        ikeep = np.where(choose < np.exp(- xpl))[0]
        if(len(ikeep) > (size - nkeep)):
            ikeep = ikeep[0:len(ikeep) - (size - nkeep)]
        xschechter[nkeep:nkeep + len(ikeep)] = xpl[ikeep]
        nkeep = nkeep + len(ikeep)
    return(xschechter)
    

class SchechterLikelihood(object):
    def __init__(self, log_lambdas=None, log_lambda_threshold=None, properties=None,
                 properties_threshold=None, model_type=None, alpha_0_range=[-2., 2.5],
                 log_lambda_min_range=[-10., -7.8], log_lambda_star_range=[-4., 2.],
                 weights=None, weights_threshold=None):
        self.log_lambdas = log_lambdas
        self.log_lambda_threshold = log_lambda_threshold
        self.properties = properties
        self.properties_threshold = properties_threshold
        if(weights is None):
            self.weights = np.ones(len(self.log_lambdas), dtype=np.float32)
        else:
            self.weights = weights
        if(weights_threshold is None):
            self.weights_threshold = np.ones(len(self.log_lambda_threshold), dtype=np.float32)
        else:
            self.weights_threshold = weights_threshold
        self.model_type = model_type
        self.alpha_0_range = alpha_0_range
        self.log_lambda_min_range = log_lambda_min_range
        self.log_lambda_star_range = log_lambda_star_range
        self.set_correction_interp()
        self.sign = 1.
        return

    def __call__(self, theta=None):
        alpha_0 = theta[0]
        log_lambda_star = theta[1]
        log_lambda_min = theta[2]
        if((alpha_0 < self.alpha_0_range[0]) |
           (alpha_0 > self.alpha_0_range[1]) |
           (log_lambda_star < self.log_lambda_star_range[0]) |
           (log_lambda_star > self.log_lambda_star_range[1]) |
           (log_lambda_min < self.log_lambda_min_range[0]) |
           (log_lambda_min > self.log_lambda_min_range[1])):
            return(self.sign * (- 1.e+30))
        if(self.model_type is not None):
            beta = theta[3]
        else:
            beta = 0.
        return(self.sign * self.loglikelihood(alpha_0=alpha_0, beta=beta,
                                              log_lambda_star=log_lambda_star,
                                              log_lambda_min=log_lambda_min))

    def alpha_prop(self, alpha_0=None, prop=None, beta=None, model_type=None):
        '''Parametrize Schechter function slope alpha as a function of some property'''
        if(model_type is None):
            prop_0 = 1.
            alpha = alpha_0 * (prop / prop_0)**beta

        return(alpha)

    def schechter_x(self, x=None, alpha=None):
        return(schechter_function(x=x, alpha=alpha))

    def power_law_int(self, log_lam_min=- 7., log_lam_max=1., alpha=None, lambda_star=None):
        """Integral of the power law component of the Schechter function"""
    
        lam_min = 10**log_lam_min
        lam_max = 10**log_lam_max

        f = (lam_max / lambda_star)**(- alpha + 1) - (lam_min / lambda_star)**(- alpha + 1)
        
        return lambda_star * f / (- alpha + 1)

    def correction_func(self, x, alpha):
        """Correction term between power law and Schechter function"""
        f = np.log(10.)*(1. - np.exp(- x)) * x**(- alpha + 1)
        return f

    def set_correction_interp(self, n=1000, nalphas=400, alpha_min=-3.,
                              alpha_max=3.5, log_x_lo=-10., log_x_hi=6.):
        """Creates 2D interpolation of the correction term integral term for log_lam and alpha grid"""

        alphas = np.linspace(alpha_min, alpha_max, nalphas)
        
        d_log_x = (log_x_hi - log_x_lo) / n
        log_x = log_x_lo + d_log_x * (np.arange(n, dtype=np.float64) + 0.5)
        log_xs = log_x + 0.5 * d_log_x

        alpalp, xx = np.meshgrid(alphas, log_xs)

        phi = self.correction_func(10**xx, alpalp)
        z = np.cumsum(phi, axis=0) * d_log_x
        
        interp_obj = scipy.interpolate.RectBivariateSpline(10**log_xs, alphas, z)

        self.correction_interp = interp_obj
        return

    def total_integral(self, log_lam_min=- 7.,log_lam_max=1., alpha=None, lambda_star=None):
        """Calculates total integral"""
        lam_min = 10**log_lam_min
        lam_max = 10**log_lam_max

        alpha_unique, alpha_index, alpha_count = np.unique(alpha, return_index=True,
                                                           return_counts=True)

        if(len(alpha_unique) == 1):
            if((- alpha_unique[0] + 1.) == 0):
                alpha_unique[0] = 1.000001
            I1 = np.zeros(len(alpha), dtype=np.float32)
            I2 = np.zeros(len(alpha), dtype=np.float32)
            I1[:] = self.power_law_int(log_lam_min, log_lam_max, alpha_unique, lambda_star)
            I2[:] = lambda_star * (self.correction_interp(lam_max / lambda_star, alpha_unique, grid=False) -
                                   self.correction_interp(lam_min / lambda_star, alpha_unique, grid=False))

            ibad = np.where(I1 != I1)[0]
            if(len(ibad) > 0):
                print(lam_max)
                print(lambda_star)
                print(alpha_unique)
                print("I1[bad] = {i}".format(i=I1[ibad]))

            ibad = np.where(I2 != I2)[0]
            if(len(ibad) > 0):
                print(lam_max)
                print(lambda_star)
                print(alpha_unique)
                print("I2[bad] = {i}".format(i=I2[ibad]))
        else:
            ione = np.where((- alpha + 1.) == 0.)[0]
            alpha[ione] = 1.000001
            I1 = self.power_law_int(log_lam_min, log_lam_max, alpha, lambda_star)
            I2 = lambda_star * (self.correction_interp(lam_max / lambda_star, alpha, grid=False) -
                                self.correction_interp(lam_min / lambda_star, alpha, grid=False))

        return I1 - I2

    def gamma_integral(self, log_min=None, log_max=None, alpha=None, lambda_star=None):
        """Integral of the Schehter function using Gamma function, value for alpha<-1"""

        alpha_unique, alpha_index, alpha_count = np.unique(alpha, return_index=True,
                                                           return_counts=True)

        if(len(alpha_unique) == 1):
            f = np.zeros(len(alpha), dtype=np.float32)
            tmpf = (lambda_star * (scipy.special.gammaincc(np.float64(1. - alpha_unique[0]), np.float64((10.**log_min) / lambda_star)) -
                                   scipy.special.gammaincc(np.float64(1. - alpha_unique[0]), np.float64((10.**log_max) / lambda_star))) *
                    scipy.special.gamma(1. - np.float64(alpha_unique[0])))
            f[:] = tmpf
        else:
            f = (lambda_star * (scipy.special.gammaincc(1. - alpha, (10.**log_min) / lambda_star) -
                                scipy.special.gammaincc(1. - alpha, (10.**log_max) / lambda_star)) *
                 scipy.special.gamma(1. - alpha))

        return f

    def log_schechter_prop(self, log_lambdas=None, prop=None, alpha_0=None, beta=None,
                           log_lambda_star=None, model_type=None):
        """Schechter as a function of properties"""
        
        alpha = self.alpha_prop(alpha_0, prop, beta, model_type)

        ln10 = np.log(10.)

        p_Xi = (np.log(ln10) + log_lambda_star * ln10 +
                (- alpha + 1.) * (log_lambdas - log_lambda_star) * ln10 -
                10.**(log_lambdas - log_lambda_star))
        return p_Xi

    def normalisation_prop(self, prop=None, alpha_0=None, beta=None,
                           log_lambda_star=None, log_lambda_min=None, model_type=None):
        """This is the normalisation constant for the schechter_mass() function."""
        integral = scipy.integrate.quad(lambda x: np.exp(self.log_schechter_prop(x, prop, alpha_0, beta,
                                                                                 log_lambda_star,
                                                                                 model_type)),
                                        log_lambda_min,log_lambda_max)[0]

        return 1. / integral

    def loglikelihood(self, alpha_0=None, beta=None, log_lambda_star=None,
                      log_lambda_min=None):

        lambda_star = 10.**log_lambda_star
        
        alphas = self.alpha_prop(alpha_0=alpha_0, prop=self.properties, beta=beta,
                                 model_type=self.model_type)
        
        inds_alp_lt_1 = np.where(alphas < 1)
        inds_alp_gt_1 = np.where(alphas >= 1)

        norms = np.ones(len(alphas))

        norms[inds_alp_gt_1] = 1. / (self.total_integral(log_lambda_min, 2.2,
                                                         alphas[inds_alp_gt_1],
                                                         lambda_star))
        norms[inds_alp_lt_1] = 1. / (self.gamma_integral(log_lambda_min, 2.2,
                                                         alphas[inds_alp_lt_1],
                                                         lambda_star))

        log_likelihoods = self.log_schechter_prop(self.log_lambdas,
                                                  self.properties,
                                                  alpha_0,
                                                  beta,
                                                  log_lambda_star,
                                                  self.model_type)

        ln_likelihood = ((log_likelihoods + np.log(norms)) * self.weights).sum()

        alphas_threshold = self.alpha_prop(alpha_0, self.properties_threshold, beta,
                                           self.model_type)
        norms_threshold = np.ones(len(alphas_threshold))

        inds_alps_lt_1 = np.where(alphas_threshold < 1)
        inds_alps_gt_1 = np.where(alphas_threshold > 1)

        norms_threshold[inds_alps_gt_1] = 1. / (self.total_integral(log_lambda_min,
                                                                    2.2,
                                                                    alphas_threshold[inds_alps_gt_1],
                                                                    lambda_star))
        norms_threshold[inds_alps_lt_1] = 1. / (self.gamma_integral(log_lambda_min,
                                                                    2.2,
                                                                    alphas_threshold[inds_alps_lt_1],
                                                                    lambda_star))

        likelihoods_alps_gt_1 = (norms_threshold[inds_alps_gt_1] *
                                 self.total_integral(log_lambda_min,
                                                     self.log_lambda_threshold[inds_alps_gt_1],
                                                     alphas_threshold[inds_alps_gt_1],
                                                     lambda_star))

        likelihoods_alps_lt_1 = (norms_threshold[inds_alps_lt_1] *
                                 self.gamma_integral(log_lambda_min,
                                                     self.log_lambda_threshold[inds_alps_lt_1],
                                                     alphas_threshold[inds_alps_lt_1],
                                                     lambda_star))

        izero = np.where(likelihoods_alps_gt_1 <= 0.)[0]
        if(len(izero) > 0):
            print("llm={l}".format(l=log_lambda_min))
            print(self.log_lambda_threshold[inds_alps_gt_1][izero], flush=True)

        izero = np.where(likelihoods_alps_lt_1 <= 0.)[0]
        if(len(izero) > 0):
            izero = np.where(likelihoods_alps_lt_1 <= 0.)[0]
            print("{n} likelihoods exactly zero".format(n=len(izero)), flush=True)
            print(likelihoods_alps_lt_1[izero], flush=True)
            print(log_lambda_min, flush=True)
            print(alpha_0, flush=True)
            print(lambda_star, flush=True)
            print(self.log_lambda_threshold[inds_alps_lt_1][izero], flush=True)
            likelihoods_alps_lt_1[izero] = 1.e-30

        ln_likelihood = (ln_likelihood +
                         (self.weights_threshold[inds_alps_gt_1] * np.log(likelihoods_alps_gt_1)).sum() +
                         (self.weights_threshold[inds_alps_lt_1] * np.log(likelihoods_alps_lt_1)).sum())

        if(np.isnan(ln_likelihood) == True):
            print("NaN!!!")
            print(alpha_0)
            print(beta)
            print(log_lambda_min)
            print(log_lambda_star)
            print(ln_likelihood)
            print(log_likelihoods.sum())
            print(log_likelihoods)
            print(norms.min())
            print(norms.max())
            if(len(likelihoods_alps_gt_1) > 0):
                print(likelihoods_alps_gt_1.min())
                print(likelihoods_alps_gt_1.max())
            if(len(likelihoods_alps_lt_1) > 0):
                print(likelihoods_alps_lt_1.min())
                print(likelihoods_alps_lt_1.max())
            print("Uh ooh", flush=True)

        if(np.isfinite(ln_likelihood) == False):
            print("Infinite!!", flush=True)
            ii = np.where(log_likelihoods <= 0)[0]
            print(self.log_lambdas[ii])
            print(self.log_lambdas)
            print(alpha_0)
            print(alphas)
            print(beta)
            print(log_lambda_min)
            print(log_lambda_star)
            print(ln_likelihood)
            print(log_likelihoods.sum())
            print(log_likelihoods)
            print(norms.min())
            print(norms.max())
            if(len(likelihoods_alps_gt_1) > 0):
                print(likelihoods_alps_gt_1.min())
                print(likelihoods_alps_gt_1.max())
            if(len(likelihoods_alps_lt_1) > 0):
                print(likelihoods_alps_lt_1.min())
                print(likelihoods_alps_lt_1.max())
            print("Uh ooh", flush=True)

        return(ln_likelihood)


def schechter_emcee(name='lambda', log_lambdas=None, log_lambda_threshold=None,
                    properties=None, properties_threshold=None,
                    weights=None, weights_threshold=None,
                    alpha_0_range=None, log_lambda_star_range=None,
                    log_lambda_min_range=None, Fagn_lambda_c_range=[-4., 0.],
                    nFagn=21, offset=0.):

    # Create likelihood object
    schlike = SchechterLikelihood(log_lambdas=log_lambdas,
                                  log_lambda_threshold=log_lambda_threshold,
                                  properties=properties,
                                  properties_threshold=properties_threshold,
                                  weights=weights,
                                  weights_threshold=weights_threshold,
                                  alpha_0_range=alpha_0_range,
                                  log_lambda_star_range=log_lambda_star_range,
                                  log_lambda_min_range=log_lambda_min_range)
    
    # Now run emcee
    nwalkers = 16
    ndim = 3
    log_lambda_star_st = - 0.5 + np.random.normal() * 0.2
    alpha_st = 1.5 + np.random.normal() * 0.2
    log_lambda_min_st = log_lambda_min_range[1] - 3. + np.random.normal() * 0.2
    theta_st = np.array([alpha_st, log_lambda_star_st, log_lambda_min_st])
    theta_0 = (np.outer(np.ones(nwalkers), theta_st) +
               np.random.normal(size=(nwalkers, ndim)) * 0.2)
    itoohigh = np.where(theta_0[:, 2] > log_lambda_min_range[1])[0]
    theta_0[itoohigh, 2] = log_lambda_min_range[1] - 0.1
    schlike.sign = 1.  # maximize

    with multiprocessing.Pool(8) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, schlike, pool=pool)

        print("emcee Burn-in", flush=True)
        state = sampler.run_mcmc(theta_0, 100)
        sampler.reset()

        print("emcee Run", flush=True)
        sampler.run_mcmc(state, 50000)

    chain_vals = sampler.get_chain(flat=True)

    print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)), flush=True)
    print("Mean autocorrelation time: {0:.3f}".format(np.mean(sampler.get_autocorr_time())), flush=True)
    
    Fagn_lambda_cs = (Fagn_lambda_c_range[0] + (Fagn_lambda_c_range[1] - Fagn_lambda_c_range[0]) *
                      np.arange(nFagn, dtype=np.float32) / np.float32(nFagn - 1))

    if(chain_vals.shape[0] < 10000):
        nchain = chain_vals.shape[0]
    else:
        nchain = 10000
    
    fit_dtype = np.dtype([('loglike', np.float32),
                          ('alpha', np.float32),
                          ('log_{n}_star'.format(n=name), np.float32),
                          ('log_{n}_min'.format(n=name), np.float32),
                          ('Fagn', np.float32, nFagn)])
    chain = np.zeros(nchain, dtype=fit_dtype)
    
    for i in np.arange(len(chain), dtype=np.int32):
        if((i % 1000) == 0):
            print(i, flush=True)
        k = chain_vals.shape[0] - 1 - i
        chain['alpha'][i] = chain_vals[k, 0]
        chain['log_{n}_star'.format(n=name)][i] = chain_vals[k, 1] + offset
        chain['log_{n}_min'.format(n=name)][i] = chain_vals[k, 2] + offset
        theta = chain_vals[k, :]
        chain['loglike'][i] = schlike(theta)
        for j, Fagn_lambda_c in enumerate(Fagn_lambda_cs):
            chain['Fagn'][i, j] = Fabove(alpha=chain['alpha'][i],
                                         log_lambda_star=chain['log_{n}_star'.format(n=name)][i] - offset,
                                         log_lambda_min=chain['log_{n}_min'.format(n=name)][i] - offset,
                                         log_lambda_c=Fagn_lambda_c - offset)

    return(chain)
