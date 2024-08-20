import numpy as np
import scipy.interpolate
import scipy.integrate
import scipy.special


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
    def __init__(self, log_lambdas=None, log_lambda_limits=None, properties=None,
                 properties_limits=None, model_type=None, alpha_0_range=[-2., 2.5],
                 log_lambda_min_range=[-10., -3.5], log_lambda_star_range=[-4., 2.]):
        self.log_lambdas = log_lambdas
        self.log_lambda_limits = log_lambda_limits
        self.properties = properties
        self.properties_limits = properties_limits
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

    def set_correction_interp(self, n=1000, nalphas=100, alpha_min=-3.,
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
            I1 = np.zeros(len(alpha), dtype=np.float32)
            I2 = np.zeros(len(alpha), dtype=np.float32)
            I1[:] = self.power_law_int(log_lam_min, log_lam_max, alpha_unique, lambda_star)
            I2[:] = lambda_star * (self.correction_interp(lam_max / lambda_star, alpha_unique, grid=False) -
                                   self.correction_interp(lam_min / lambda_star, alpha_unique, grid=False))
        else:
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
            f[:] = (lambda_star * (scipy.special.gammainc(1. - alpha_unique, (10.**log_max) / lambda_star) -
                                   scipy.special.gammainc(1. - alpha_unique, (10.**log_min) / lambda_star)) *
                    scipy.special.gamma(1. - alpha_unique))
        else:
            f = (lambda_star * (scipy.special.gammainc(1. - alpha, (10.**log_max) / lambda_star) -
                                scipy.special.gammainc(1. - alpha, (10.**log_min) / lambda_star)) *
                 scipy.special.gamma(1. - alpha))

        return f

    def schechter_prop(self, lambdas=None, prop=None, alpha_0=None, beta=None,
                       lambda_star=None, model_type=None):
        """Schechter as a function of properties"""
        
        alpha = self.alpha_prop(alpha_0, prop, beta, model_type)

        p_Xi = (np.log(10.) * lambda_star * (lambdas / lambda_star)**(- alpha + 1.) *
                np.exp( - lambdas / lambda_star))
        return p_Xi

    def normalisation_prop(self, prop=None, alpha_0=None, beta=None,
                           lambda_star=None, log_lambda_min=None, model_type=None):
        """This is the normalisation constant for the schechter_mass() function."""
        integral = scipy.integrate.quad(lambda x: schechter_prop(10.**x, prop, alpha_0, beta,
                                                                 lambda_star, model_type),
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

        likelihoods = norms * self.schechter_prop(10.**self.log_lambdas,
                                                  self.properties,
                                                  alpha_0,
                                                  beta,
                                                  lambda_star,
                                                  self.model_type)

        ln_likelihood = (np.log(likelihoods)).sum()

        alphas_limits = self.alpha_prop(alpha_0, self.properties_limits, beta,
                                        self.model_type)
        norms_limits = np.ones(len(alphas_limits))

        inds_alps_lt_1 = np.where(alphas_limits < 1)
        inds_alps_gt_1 = np.where(alphas_limits > 1)

        norms_limits[inds_alps_gt_1] = 1. / (self.total_integral(log_lambda_min,
                                                                 2.2,
                                                                 alphas_limits[inds_alps_gt_1],
                                                                 lambda_star))
        norms_limits[inds_alps_lt_1] = 1. / (self.gamma_integral(log_lambda_min,
                                                                 2.2,
                                                                 alphas_limits[inds_alps_lt_1],
                                                                 lambda_star))

        likelihoods_alps_gt_1 = (norms_limits[inds_alps_gt_1] *
                                 self.total_integral(log_lambda_min,
                                                     self.log_lambda_limits[inds_alps_gt_1],
                                                     alphas_limits[inds_alps_gt_1],
                                                     lambda_star))

        likelihoods_alps_lt_1 = (norms_limits[inds_alps_lt_1] *
                                 self.gamma_integral(log_lambda_min,
                                                     self.log_lambda_limits[inds_alps_lt_1],
                                                     alphas_limits[inds_alps_lt_1],
                                                     lambda_star))

        ln_likelihood = (ln_likelihood +
                         (np.log(likelihoods_alps_gt_1)).sum() +
                         (np.log(likelihoods_alps_lt_1)).sum())

        return(ln_likelihood)
