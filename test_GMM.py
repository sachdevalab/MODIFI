import numpy as np
import warnings
from sklearn.mixture import GaussianMixture
from sklearn.exceptions import ConvergenceWarning


# Fit a gaussian mixture model to the ipd of a single molecule, return the binarized values
def fitGaussian(ipds, initmean = False, zmw=-1):
    with np.errstate(invalid='ignore'):
        ridx = np.where(ipds > 0)[0]

    # lrnn = np.log(ipds[ridx]).reshape(-1,1)
    lrnn = ipds[ridx].reshape(-1,1)
    print ("len(lrnn)", len(lrnn))
    gmm2 = GaussianMixture(2, covariance_type='spherical')
    gmm2.fit(lrnn)
    print ("xxx", gmm2.means_, gmm2.weights_, gmm2.converged_)

    gmm1 = GaussianMixture(1, covariance_type='spherical')
    if initmean:
        gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                              weights_init=np.array([.85, .15]), tol=1e-6)
    else:
        gmm2 = GaussianMixture(2, covariance_type='spherical')
    
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        gmm1.fit(lrnn)
        if not gmm1.converged_:
            print('zmw #%d did not converge on gmm1' % (zmw))
        gmm2.fit(lrnn)
        print (gmm2.means_, gmm2.weights_, gmm2.converged_)
        convround = 0
        if not gmm2.converged_:
            gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                  n_init=10, init_params='random', tol=1e-5, max_iter=200)
            gmm2.fit(lrnn)
            convround = 1
            if not gmm2.converged_:
                gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                      n_init=20, init_params='random', tol=1e-4, max_iter=400)
                gmm2.fit(lrnn)
                convround = 2
                if not gmm2.converged_:
                    gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                      n_init=40, init_params='random', tol=1e-3, max_iter=600)
                    gmm2.fit(lrnn)
                    convround = 3
                    if not gmm2.converged_:
                        convround = 9
                        print('zmw #%d did not converge on gmm2 even after extensions' % (zmw))
    print (gmm1.aic(lrnn), gmm2.aic(lrnn))
    aicdif = gmm1.aic(lrnn) - gmm2.aic(lrnn)
    mixmeans = gmm2.means_.flatten()
    mixweights = gmm2.weights_.flatten()
    elevstate = np.argmax(mixmeans)
    resp = gmm2.predict_proba(np.log(ipds[ridx].reshape(-1,1)))
    respfull = np.empty(len(ipds), dtype='float32')
    respfull.fill(np.nan)
    respfull[ridx] = resp[:,elevstate]
    convergInf = str(convround) + '.' + str(gmm2.n_iter_)
    print ("////")
    print (mixmeans, mixweights, aicdif, convergInf)
    ## print whether the model is 1 GMM or 2 GMM
    if aicdif >= 0:
        print ('2 GMM')
    else:
        print('1 GMM')
    
    return (respfull, np.array([mixmeans[1-elevstate], mixmeans[elevstate]]), 
            np.array([mixweights[1-elevstate], mixweights[elevstate]]), aicdif, convergInf)


# ## test the function
file = "/home/shuaiw/borg/bench/C227/WGA2/ipd_ratio/CP011331.1.ipd3.csv"
# file = "/home/shuaiw/borg/bench/breakdb/ipd_ratio/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_494_C_0_547023.ipd3.csv"
## read the ipd file using pd
import pandas as pd
df = pd.read_csv(file, nrows=10000)
X2 = df['ipd_ratio'].values.reshape(-1,1)


# ## generate random data with one or two gaussians
# np.random.seed(0)
# n_samples = 10000
# # generate random data with one gaussian
# # Generate samples from a normal distribution
# X1 = np.random.randn(n_samples) * 0.5 + 1
# X1 = X1.reshape(-1,1)

# X2 = np.random.randn(100) * 1 + 5
# X2 = X2.reshape(-1,1)

# # # generate random data with two gaussians
# X2 = np.vstack([X1, X2])
# print (X2)
# # ## fit the model
# respfull, mixmeans, mixweights, aicdif, convergInf = fitGaussian(X2.flatten(), True)
# # print (respfull, mixmeans, mixweights, aicdif, convergInf)
# print (aicdif)
gmm2 = GaussianMixture(2, weights_init=np.array([.99, .01]), means_init=np.array([1.1, 3]).reshape((2,1)))
gmm2.fit(X2)
print (gmm2.means_, gmm2.aic(X2), gmm2.weights_.flatten())

gmm1 = GaussianMixture(1)
gmm1.fit(X2)
print (gmm1.means_, gmm1.aic(X2))
## plot the distribution of X2
import matplotlib.pyplot as plt
import seaborn as sns
sns.histplot(X2)
## save the plot
plt.savefig('tmp/test.png')