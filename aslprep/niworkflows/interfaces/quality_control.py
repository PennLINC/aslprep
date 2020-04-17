import numpy as np
import nibabel as nib
import pandas as pd 
from nibabel.processing import smooth_image
from scipy.stats import gmean

def dc(input1, input2):
    r"""
    Dice coefficient
    Computes the Dice coefficient (also known as Sorensen index) between the binary
    objects in two images.
    The metric is defined as
    .. math::
        DC=\frac{2|A\cap B|}{|A|+|B|}
    , where :math:`A` is the first and :math:`B` the second set of samples (here: binary objects).
    Parameters
    ----------
    input1 : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    input2 : array_like
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    Returns
    -------
    dc : float
        The Dice coefficient between the object(s) in ```input1``` and the
        object(s) in ```input2```. It ranges from 0 (no overlap) to 1 (perfect overlap).
    Notes
    -----
    This is a real metric.
    """
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))

    intersection = np.count_nonzero(input1 & input2)

    size_i1 = np.count_nonzero(input1)
    size_i2 = np.count_nonzero(input2)

    try:
        dc = 2. * intersection / float(size_i1 + size_i2)
    except ZeroDivisionError:
        dc = 0.0

    return dc


def jc(input1, input2):
    r"""
    Jaccard coefficient
    Computes the Jaccard coefficient between the binary objects in two images.
    Parameters
    ----------
    input1: array_like
            Input data containing objects. Can be any type but will be converted
            into binary: background where 0, object everywhere else.
    input2: array_like
            Input data containing objects. Can be any type but will be converted
            into binary: background where 0, object everywhere else.
    Returns
    -------
    jc: float
        The Jaccard coefficient between the object(s) in `input1` and the
        object(s) in `input2`. It ranges from 0 (no overlap) to 1 (perfect overlap).
    Notes
    -----
    This is a real metric.
    """
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))

    intersection = np.count_nonzero(input1 & input2)
    union = np.count_nonzero(input1 | input2)

    jc = float(intersection) / float(union)

    return jc

def crosscorr(input1,input2):

   """ 
   cross correlation
   computer compute cross correction bewteen input mask 
   """

   input1 = np.atleast_1d(input1.astype(np.bool))
   input2 = np.atleast_1d(input2.astype(np.bool))

   from scipy.stats.stats import pearsonr 
   cc=pearsonr(input1,input2)
   return cc 

def coverage(input1,input2):
    """
    estimate the coverage between  two mask
    """
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool)) 
     
    intsec=np.count_nonzero(input1 & input2)
    if np.sum(input1)> np.sum(input2):
        smallv=np.sum(input2)
    else:
        smallv=np.sum(input1)
    cov=float(intsec)/float(smallv)
    return cov


def cbf_qei(gm,wm,csf,img,thresh=0.7):
    
    def fun1(x,xdata):
        d1=np.exp(-(x[0])*np.power(xdata,x[1]))
        return(d1)
    
    def fun2(x,xdata):
        d1=1-np.exp(-(x[0])*np.power(xdata,x[1]))
        return(d1)
    
    x1 = [0.054,0.9272]; x2 = [2.8478,0.5196]; x4 = [3.0126, 2.4419]                                     
    scbf=smooth_image(nib.load(img),fwhm=5).get_fdata()# smooth the image 
    #load prob maps
    gmm=nib.load(gm).get_fdata(); wmm=nib.load(wm).get_fdata(); ccf=nib.load(csf).get_fdata()
    pbcf=2.5*gmm+wmm  # gmm is 2.5 times wm
    msk=np.array((scbf!= 0)&(scbf != np.nan )&(pbcf != np.nan )).astype(int)

    gm1=np.array(gmm>thresh)
    wm1=np.array(wmm>thresh)
    cc1=np.array(ccf>thresh)
    r1=np.array([0,np.corrcoef(scbf[msk==1],pbcf[msk==1])[1,0]]).max()
   
    V=((np.sum(gm1)-1)*np.var(scbf[gm1>0])+(np.sum(wm1)-1)*np.var(scbf[wm1>0])+(np.sum(cc1)-1)  \
           *np.var(scbf[cc1>0]))/(np.sum(gm1>0)+np.sum(wm1>0)+np.sum(cc1>0)-3)
    
    negGM=np.sum(scbf[gm1]<0)/(np.sum(gm1))
    GMCBF=np.mean(scbf[gm1])
    CV=V/np.abs(GMCBF)
    Q = [fun1(x1,CV),fun1(x2,negGM),fun2(x4,r1)]
    return gmean(Q)


