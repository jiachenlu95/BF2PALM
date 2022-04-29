import numpy as np
def cal_variance(domain):
    profile = np.zeros(domain.shape[0])
    for i in range(domain.shape[0]):
        profile[i] = np.nanmean(domain[i,:])
        
    return(profile,np.var(profile))