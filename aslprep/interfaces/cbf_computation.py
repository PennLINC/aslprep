import os
import numpy as np
import pandas as pd
import nibabel as nb
from nibabel.processing import smooth_image
from scipy.stats import gmean
from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import (traits, TraitedSpec, BaseInterfaceInputSpec, SimpleInterface,
                                    File, isdefined)
from nipype.interfaces.fsl.base import (FSLCommand, FSLCommandInputSpec)
from nipype.interfaces.ants import ApplyTransforms
from pkg_resources import resource_filename as pkgrf

LOGGER = logging.getLogger('nipype.interface')


class _refinemaskInputSpec(BaseInterfaceInputSpec):
    in_t1mask = File(exists=True, mandatory=True, desc='t1 mask')
    in_aslmask = File(exists=True, mandatory=True, desct='asl mask')
    transforms = File(exists=True, mandatory=True, desc='transfom')
    out_mask = File(exists=False, mandatory=False, desc='output mask')
    out_tmp = File(exists=False, mandatory=False, desc='tmp mask')


class _refinemaskOutputSpec(TraitedSpec):
    out_mask = File(exists=False, desc='output mask')
    out_tmp = File(exists=False, desc='tmp mask')


class refinemask(SimpleInterface):
    r"""
    the code refine the asl mask with t1w mask
    the output is refined asl mask

    """
    input_spec = _refinemaskInputSpec
    output_spec = _refinemaskOutputSpec

    def _run_interface(self, runtime):
        self._results['out_tmp'] = fname_presuffix(self.inputs.in_aslmask,
                                                   suffix='_tempmask', newpath=runtime.cwd)
        self._results['out_mask'] = fname_presuffix(self.inputs.in_aslmask,
                                                    suffix='_refinemask', newpath=runtime.cwd)
        b1 = ApplyTransforms()
        b1.inputs.dimension = 3
        b1.inputs.float = True
        b1.inputs.input_image = self.inputs.in_t1mask
        b1.inputs.interpolation = 'NearestNeighbor'
        b1.inputs.reference_image = self.inputs.in_aslmask
        b1.inputs.transforms = self.inputs.transforms
        b1.inputs.input_image_type = 3
        b1.inputs.output_image = self._results['out_tmp']
        b1.run()

        from nipype.interfaces.fsl import MultiImageMaths
        mat1 = MultiImageMaths()
        mat1.inputs.in_file = self._results['out_tmp']
        mat1.inputs.op_string = " -mul  %s -bin"
        mat1.inputs.operand_files = self.inputs.in_aslmask
        mat1.inputs.out_file = self._results['out_mask']
        mat1.run()
        self.inputs.out_mask = os.path.abspath(self._results['out_mask'])
        return runtime


class _extractCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='raw asl file')
    asl_file = File(exists=True, mandatory=True, desc='preprocessed asl file')
    in_mask = File(exists=True, mandatory=True, desc='mask')
    dummy_vols = traits.Int(default_value=0, exit=False, mandatory=False,
                            desc='remove first n volumes')
    in_metadata = traits.Dict(exists=True, mandatory=True,
                              desc='metadata for asl or deltam ')
    bids_dir=traits.Str(exits=True,mandatory=True,desc=' bids directory')
    fwhm = traits.Float(default_value=5, exists=True, mandatory=False, desc='fwhm')
    out_file = File(exists=False, mandatory=False, desc='cbf timeries data')
    out_avg = File(exists=False, mandatory=False, desc='average control')


class _extractCBFOutputSpec(TraitedSpec):
    out_file = File(exists=False, desc='cbf timeries data')
    out_avg = File(exists=False, desc='average control')


class extractCBF(SimpleInterface):
    """
    extract  CBF timeseries
    by substracting label from control
    or viceversa

    it generate M0 maps for cbf computation

    """

    input_spec = _extractCBFInputSpec
    output_spec = _extractCBFOutputSpec

    def _run_interface(self, runtime):
        file1 = os.path.abspath(self.inputs.in_file)
        # check if there is m0 file
        m0num=1
        m0file=[]
        aslfile_linkedM0=[]

        if self.inputs.in_metadata['M0'] != "True" and self.inputs.in_metadata['M0'] != "False" and type(self.inputs.in_metadata['M0']) != int :
            m0file=os.path.abspath(self.inputs.bids_dir+'/'+self.inputs.in_metadata['M0'])
            m0file_metadata=readjson(m0file.replace('nii.gz','json'))
            aslfile_linkedM0 = os.path.abspath(self.inputs.bids_dir+'/'+m0file_metadata['IntendedFor'])
        elif type(self.inputs.in_metadata['M0']) == int :
            m0num=int(self.inputs.in_metadata['M0'])
        else:
            print('no M0 file or numerical M0, the average control will be used \
             in the case of deltam, M0 is required for cbf quantifcation') 
        
        aslcontext1 = file1.replace('_asl.nii.gz', '_aslcontext.tsv')
        aslcontext = pd.read_csv(aslcontext1)
        idasl = aslcontext['volume_type'].tolist()

        # get the control,tag,moscan or label 
        controllist = [i for i in range(0, len(idasl)) if idasl[i] == 'control']
        labellist = [i for i in range(0, len(idasl)) if idasl[i] == 'label']
        m0list = [i for i in range(0, len(idasl)) if idasl[i] == 'm0scan']

        deltamlist = [i for i in range(0, len(idasl)) if idasl[i] == 'deltam']

        
        allasl = nb.load(self.inputs.asl_file)
        mask = nb.load(self.inputs.in_mask).get_fdata()
        dataasl = allasl.get_fdata()

        if len(dataasl.shape) == 5:
            raise RuntimeError('Input image (%s) is 5D.')

        if len(deltamlist) > 0 : 
            cbf_data = dataasl[:, :, :, deltamlist]
        elif len(labellist) > 0 :
            control_img = dataasl[:, :, :, controllist]
            label_img = dataasl[:, :, :, labellist] 
            cbf_data = np.subtract(control_img, label_img)
        else: 
            raise RuntimeError('no valid asl or cbf image.')
      
        
        if self.inputs.dummy_vols != 0:
            cbf_data = np.delete(cbf_data, range(0, self.inputs.dummy_vols), axis=3)
            #control_img = np.delete(control_img, range(0, self.inputs.dummy_vols), axis=3)

        # MO file
        if m0file or aslfile_linkedM0 :
            # get the raw m0 file also check intended for
            #m0file=nb.load(m0file).get_fdata()
            #regsiter m0file to aslfile here
            if m0file:
                m0file = m0file
            else:
                m0file = aslfile_linkedM0

            newm0 = fname_presuffix(self.inputs.asl_file,
                                                    suffix='_m0file') 
            newm0 = regmotoasl(asl=self.inputs.asl_file,m0file=m0file,m02asl=newm0)
            m0data_smooth = smooth_image(nb.load(newm0), fwhm=self.inputs.fwhm).get_data()
            if len(m0data_smooth.shape) > 3 :
                avg_control = mask*np.mean(m0data_smooth, axis=3)
            else:
                avg_control = mask*m0data_smooth

        elif len(m0list) > 0 and self.inputs.in_metadata['M0'] == "True" :
            # if no m0file, check from asl data
            modata2 = dataasl[:, :, :, m0list]
            con2 = nb.Nifti1Image(modata2, allasl.affine, allasl.header)
            m0data_smooth = smooth_image(con2, fwhm=self.inputs.fwhm).get_data()
            if len(m0data_smooth.shape) > 3 :
                avg_control = mask*np.mean(m0data_smooth, axis=3)
            else:
                avg_control = mask*m0data_smooth
        elif len(controllist) > 0:
            # else use average control
            control_img = dataasl[:, :, :, controllist]
            con = nb.Nifti1Image(control_img, allasl.affine, allasl.header)
            control_img1 = smooth_image(con, fwhm=self.inputs.fwhm).get_data()
            avg_control = mask*np.mean(control_img1, axis=3)
        else: 
            'precomputed m0 number will be used'
            avg_control = mask*(np.mean(np.ones_like(cbf_data),axis=3))
            avg_control = m0num*avg_control


        self._results['out_file'] = fname_presuffix(self.inputs.in_file,
                                                    suffix='_cbftimeseries', newpath=runtime.cwd)
        self._results['out_avg'] = fname_presuffix(self.inputs.in_file,
                                                   suffix='_avg_control', newpath=runtime.cwd)
        nb.Nifti1Image(
            np.divide(cbf_data,m0num), allasl.affine, allasl.header).to_filename(
            self._results['out_file'])
        nb.Nifti1Image(
            avg_control,allasl.affine, allasl.header).to_filename(
            self._results['out_avg'])

        self.inputs.out_file = os.path.abspath(self._results['out_file'])
        self.inputs.out_avg = os.path.abspath(self._results['out_avg'])
        return runtime


class _computeCBFInputSpec(BaseInterfaceInputSpec):
    in_cbf = File(exists=True, mandatory=True, desc='cbf nifti')
    in_metadata = traits.Dict(exists=True, mandatory=True,
                              desc='metadata for CBF ')
    in_m0scale=traits.Float(exists=True, mandatory=True,
                              desc='relative scale between asl and m0')
    in_m0file = File(exists=True, mandatory=False, desc='M0 nifti file')
    in_mask = File(exists=True, mandatory=False, desc='mask')
    out_cbf = File(exists=False, mandatory=False, desc='cbf timeries data')
    out_mean = File(exists=False, mandatory=False, desc='average control')
    out_att = File(exists=False, mandatory=False, desc='Arterial Transit Time')


class _computeCBFOutputSpec(TraitedSpec):
    out_cbf = File(exists=False, desc='cbf timeries data')
    out_mean = File(exists=False, desc='average control')
    out_att = File(exists=False, desc='Arterial Transit Time')


class computeCBF(SimpleInterface):
    """
    compute cbf pASL or pCASL
    """
    input_spec = _computeCBFInputSpec
    output_spec = _computeCBFOutputSpec

    def _run_interface(self, runtime):
        cbf, meancbf, att = cbfcomputation(metadata=self.inputs.in_metadata,m0scale=self.inputs.in_m0scale,
                                           mask=self.inputs.in_mask, m0file=self.inputs.in_m0file,
                                           cbffile=self.inputs.in_cbf)
        self._results['out_cbf'] = fname_presuffix(self.inputs.in_cbf,
                                                   suffix='_cbf', newpath=runtime.cwd)
        self._results['out_mean'] = fname_presuffix(self.inputs.in_cbf,
                                                    suffix='_meancbf', newpath=runtime.cwd)
        samplecbf = nb.load(self.inputs.in_m0file)
        nb.Nifti1Image(
            cbf, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_cbf'])
        nb.Nifti1Image(
            meancbf, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_mean'])
        if att is not None:
            self._results['out_att'] = fname_presuffix(self.inputs.in_cbf,
                                                       suffix='_att', newpath=runtime.cwd)
            nb.Nifti1Image(att, samplecbf.affine, samplecbf.header).to_filename(
                           self._results['out_att'])
            self.inputs.out_att = os.path.abspath(self._results['out_att'])
        self.inputs.out_cbf = os.path.abspath(self._results['out_cbf'])
        self.inputs.out_mean = os.path.abspath(self._results['out_mean'])
        # we dont know why not zeros background $
        from nipype.interfaces.fsl import MultiImageMaths
        mat1 = MultiImageMaths()
        mat1.inputs.in_file = self.inputs.out_mean
        mat1.inputs.op_string = " -mul  %s "
        mat1.inputs.operand_files = self.inputs.in_mask
        mat1.inputs.out_file = self.inputs.out_mean
        mat1.run()
        mat1 = MultiImageMaths()
        mat1.inputs.in_file = self.inputs.out_cbf
        mat1.inputs.op_string = " -mul  %s "
        mat1.inputs.operand_files = self.inputs.in_mask
        mat1.inputs.out_file = self.inputs.out_cbf
        mat1.run()
        return runtime

def cbfcomputation(metadata, mask, m0file, cbffile, m0scale=1):
    
    """
    compute cbf with pld and multi pld
    metadata
      cbf metadata
    mask
      asl mask in native space
    m0file
      m0scan
    cbffile
      already processed cbf  after tag-control substraction
    m0scale
      relative scale between m0scan and asl, default is 1
    """
    labeltype = metadata['LabelingType']
    tau = metadata['LabelingDuration']
    plds = np.array(metadata['PostLabelingDelay'])
    #m0scale = metadata['M0']
    magstrength = metadata['MagneticFieldStrength']
    t1blood = (110*int(magstrength)+1316)/1000 # https://onlinelibrary.wiley.com/doi/pdf/10.1002/mrm.24550 
    inverstiontime = np.add(tau, plds)
    #mask = nb.load(mask).get_fdata()

        
    if 'LabelingEfficiency' in metadata.keys():
        labeleff = metadata['LabelingEfficiency']
    elif 'CASL' in labeltype:
        labeleff = 0.72
    elif 'PASL' in labeltype:
        labeleff = 0.8
    else:
        print('no labelelling effiecieny')
    part_coeff = 0.9   # brain partition coefficient



    if 'CASL' in labeltype:
        pf1 = (6000*part_coeff)/(2*labeleff*t1blood*(1-np.exp(-(tau/t1blood))))
        perfusion_factor = pf1*np.exp(plds/t1blood)
    elif 'PASL' in labeltype:
        pf1 = (6000*part_coeff)/(2*labeleff)
        perfusion_factor = (pf1*np.exp(inverstiontime/t1blood))/inverstiontime
    #perfusion_factor = np.array(perfusion_factor)
    #print(perfusion_factor)

    maskx = nb.load(mask).get_fdata()
    m0data = nb.load(m0file).get_fdata()
    m0data=m0data[maskx==1]
    # compute cbf
    cbf_data = nb.load(cbffile).get_fdata()
    cbf_data=cbf_data[maskx==1]
    cbf1 = np.zeros(cbf_data.shape)
    if len(cbf_data.shape) < 2: 
        cbf1 = np.divide(cbf_data,(m0scale*m0data))
    else: 
        for i in range(cbf1.shape[1]):
            cbf1[:, i] = np.divide(cbf_data[:,i], (m0scale*m0data))
        # m1=m0scale*m0_data
        # cbf1=np.divide(cbf_data,m1)
        # for compute cbf for each PLD and TI
    att = None  
    if hasattr(perfusion_factor, '__len__'):
        permfactor = np.tile(perfusion_factor ,int(cbf_data.shape[1]/len(perfusion_factor)))
        cbf_data_ts = np.zeros(cbf_data.shape)

        #calculate  cbf with multiple plds 
        for i in range(cbf_data.shape[1]):
            cbf_data_ts[:, i] =np.multiply(cbf1[:, i],permfactor[i])
        cbf = np.zeros([cbf_data_ts.shape[0], int(cbf_data.shape[1]/len(perfusion_factor))])
        cbf_xx=np.split(cbf_data_ts,int(cbf_data_ts.shape[1]/len(perfusion_factor)),axis=1)
        
        # calculate weighted cbf with multiplds
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791289/
        # https://pubmed.ncbi.nlm.nih.gov/22084006/
        for k in range(len(cbf_xx)):
            cbf_plds = cbf_xx[k]
            pldx = np.zeros([cbf_plds.shape[0],len(cbf_plds)])
            for j in range(cbf_plds.shape[1]):
                pldx[:,j] = np.array(np.multiply(cbf_plds[:,j],plds[j]))
            cbf[:, k]=np.divide(np.sum(pldx,axis=1),np.sum(plds))
    else:
        cbf = cbf1*np.array(perfusion_factor)
        # cbf is timeseries
    # return cbf to nifti shape
    if len(cbf.shape) < 2:
        tcbf=np.zeros(maskx.shape)
        tcbf[maskx==1]=cbf
    else:
        tcbf=np.zeros([maskx.shape[0],maskx.shape[1],maskx.shape[2],cbf.shape[1]])
        for i in range(cbf.shape[1]):
            tcbfx=np.zeros(maskx.shape); 
            tcbfx[maskx==1]=cbf[:,i]
            tcbf[:,:,:,i]=tcbfx
    if len(tcbf.shape) < 4:
        meancbf = tcbf
    else:
        meancbf = np.nanmean(tcbf, axis=3)
    meancbf = np.nan_to_num(meancbf)
    tcbf = np.nan_to_num(tcbf)
    att = np.nan_to_num(att)
    return tcbf, meancbf, att

# score and scrub


class _scorescrubCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True,
                   desc='computed CBF from computeCBF')
    in_greyM = File(exists=True, mandatory=True, desc='grey  matter')
    in_whiteM = File(exists=True, mandatory=True, desc='white  matter')
    in_mask = File(exists=True, mandatory=True, desc='mask')
    in_csf = File(exists=True, mandatory=True, desc='csf')
    in_thresh = traits.Float(default_value=0.7, exists=True,
                             mandatory=False, desc='threshold of propbaility matter')
    in_wfun = traits.Str(exists=True, mandatory=False, default_value='huber',
                         option=['bisquare', 'andrews', 'cauchy', 'fair', 'logistics',
                                 'ols', 'talwar', 'welsch'],
                         desc='wavelet fun ')
    out_score = File(exists=False, desc='score timeseries data')
    out_avgscore = File(exists=False, desc='average score')
    out_scrub = File(exists=False, desc='average scrub')
    out_scoreindex = File(exists=False, desc='index of volume remove or leave by score')


class _scorescrubCBFOutputSpec(TraitedSpec):
    out_score = File(exists=False, mandatory=False, desc='score timeseries data')
    out_avgscore = File(exists=False, mandatory=False, desc='average score')
    out_scrub = File(exists=False, mandatory=False, desc='average scrub')
    out_scoreindex = File(exists=False, mandatory=False, desc='index of volume remove ')


class scorescrubCBF(SimpleInterface):
    """
    compute score and scrub
    """
    input_spec = _scorescrubCBFInputSpec
    output_spec = _scorescrubCBFOutputSpec

    def _run_interface(self, runtime):
        cbf_ts = nb.load(self.inputs.in_file).get_fdata()
        mask = nb.load(self.inputs.in_mask).get_fdata()
        greym = nb.load(self.inputs.in_greyM).get_fdata()
        whitem = nb.load(self.inputs.in_whiteM).get_fdata()
        csf = nb.load(self.inputs.in_csf).get_fdata()
        if len(cbf_ts.shape) > 3:
            cbf_scorets, index_score = _getcbfscore(cbfts=cbf_ts, wm=whitem,
                                                gm=greym, csf=csf, mask=mask,
                                                thresh=self.inputs.in_thresh)
            cbfscrub = _scrubcbf(cbf_ts=cbf_scorets, gm=greym, wm=whitem, csf=csf,
                             mask=mask, wfun=self.inputs.in_wfun, thresh=self.inputs.in_thresh)
            avgscore = np.mean(cbf_scorets, axis=3)
        else:
            cbf_scorets = cbf_ts
            index_score = np.array([0])
            cbfscrub = cbf_ts
            avgscore = cbf_ts
        
        self._results['out_score'] = fname_presuffix(self.inputs.in_file,
                                                     suffix='_cbfscorets', newpath=runtime.cwd)
        self._results['out_avgscore'] = fname_presuffix(self.inputs.in_file,
                                                        suffix='_meancbfscore',
                                                        newpath=runtime.cwd)
        self._results['out_scrub'] = fname_presuffix(self.inputs.in_file,
                                                     suffix='_cbfscrub', newpath=runtime.cwd)
        self._results['out_scoreindex'] = fname_presuffix(self.inputs.in_file,
                                                          suffix='_scoreindex.txt',
                                                          newpath=runtime.cwd, use_ext=False)
        samplecbf = nb.load(self.inputs.in_mask)

        nb.Nifti1Image(dataobj=cbf_scorets, affine=samplecbf.affine, header=samplecbf.header).to_filename(self._results['out_score'])
        nb.Nifti1Image(
             dataobj=avgscore, affine=samplecbf.affine, header=samplecbf.header).to_filename(
                self._results['out_avgscore'])
        nb.Nifti1Image(
            dataobj=cbfscrub, affine=samplecbf.affine, header=samplecbf.header).to_filename(
            self._results['out_scrub'])

        np.savetxt(self._results['out_scoreindex'],index_score, delimiter=',')

        self.inputs.out_score = os.path.abspath(self._results['out_score'])
        self.inputs.out_avgscore = os.path.abspath(self._results['out_avgscore'])
        self.inputs.out_scrub = os.path.abspath(self._results['out_scrub'])
        self.inputs.out_scoreindex = os.path.abspath(self._results['out_scoreindex'])
        return runtime


def _weightfun(x, wfun='huber'):
    """"
    get weight fun and tuner

    """
    if wfun == 'andrews':
        tuner = 1.339
        weight = (np.abs(x) < np.pi)*np.sin(x)
    elif wfun == 'bisquare':
        tuner = 4.685
        weight = (np.abs(x) < 1)*np.power((1-np.power(x, 2)), 2)
    elif wfun == 'cauchy':
        tuner = 2.385
        weight = 1/(1+np.power(x, 2))
    elif wfun == 'logistic':
        tuner = 1.205
        weight == np.tanh(x)/x
    elif wfun == 'ols':
        tuner = 1
        weight = np.repeat(1, len(x))
    elif wfun == 'talwar':
        tuner = 2.795
        weight = 1*(np.abs(x) < 1)
    elif wfun == 'welsch':
        tuner = 2.985
        weight = np.exp(-(np.power(x, 2)))
    else:
        tuner = 1.345
        weight = 1/np.abs(x)
    return weight, tuner


def _tune(wfun='huber'):
    """"
    get weight fun and tuner

    """
    if wfun == 'andrews':
        tuner = 1.339
    elif wfun == 'bisquare':
        tuner = 4.685
    elif wfun == 'cauchy':
        tuner = 2.385
    elif wfun == 'logistic':
        tuner = 1.205
    elif wfun == 'ols':
        tuner = 1
    elif wfun == 'talwar':
        tuner = 2.795
    elif wfun == 'welsch':
        tuner = 2.985
    else:
        tuner = 1.345
    return tuner


def _getchisquare(n):
    a = [0.000000, 15.484663, 8.886835, 7.224733, 5.901333, 5.126189, 4.683238, 4.272937, 4.079918,
         3.731612, 3.515615, 3.459711, 3.280471, 3.078046, 3.037280, 2.990761, 2.837119, 2.795526,
         2.785189, 2.649955, 2.637642, 2.532700, 2.505253, 2.469810, 2.496135, 2.342210, 2.384975,
         2.275019, 2.244482, 2.249109, 2.271968, 2.210340, 2.179537, 2.133762, 2.174928, 2.150072,
         2.142526, 2.071512, 2.091061, 2.039329, 2.053183, 2.066396, 1.998564, 1.993568, 1.991905,
         1.981837, 1.950225, 1.938580, 1.937753, 1.882911, 1.892665, 1.960767, 1.915530, 1.847124,
         1.947374, 1.872383, 1.852023, 1.861169, 1.843109, 1.823870, 1.809643, 1.815038, 1.848064,
         1.791687, 1.768343, 1.778231, 1.779046, 1.759597, 1.774383, 1.774876, 1.751232, 1.755293,
         1.757028, 1.751388, 1.739384, 1.716395, 1.730631, 1.718389, 1.693839, 1.696862, 1.691245,
         1.682541, 1.702515, 1.700991, 1.674607, 1.669986, 1.688864, 1.653713, 1.641309, 1.648462,
         1.630380, 1.634156, 1.660821, 1.625298, 1.643779, 1.631554, 1.643987, 1.624604, 1.606314,
         1.609462]
    b = [0, 2.177715, 1.446966, 1.272340, 1.190646, 1.151953, 1.122953, 1.103451, 1.089395,
         1.079783, 1.071751, 1.063096, 1.058524, 1.054137, 1.049783, 1.046265, 1.043192,
         1.039536, 1.038500, 1.037296, 1.033765, 1.032317, 1.031334, 1.029551, 1.028829,
         1.027734, 1.024896, 1.024860, 1.025207, 1.024154, 1.022032, 1.021962, 1.021514,
         1.020388, 1.019238, 1.020381, 1.019068, 1.018729, 1.018395, 1.017134, 1.016539,
         1.015676, 1.015641, 1.015398, 1.015481, 1.015566, 1.014620, 1.014342, 1.013901,
         1.013867, 1.013838, 1.013602, 1.013322, 1.012083, 1.013168, 1.012667, 1.011087,
         1.011959, 1.011670, 1.011494, 1.010463, 1.010269, 1.010393, 1.010004, 1.010775,
         1.009399, 1.011000, 1.010364, 1.009831, 1.009563, 1.010085, 1.009149, 1.008444,
         1.009455, 1.009705, 1.008597, 1.008644, 1.008051, 1.008085, 1.008550, 1.008265,
         1.009141, 1.008235, 1.008002, 1.008007, 1.007660, 1.007993, 1.007184, 1.008093,
         1.007816, 1.007770, 1.007932, 1.007819, 1.007063, 1.006712, 1.006752, 1.006703,
         1.006650, 1.006743, 1.007087]
    return a[n-1], b[n-1]


def _getcbfscore(cbfts, wm, gm, csf, mask, thresh=0.7):
    """ 
    score algorithm by Sudipto
    removing noisy cbf volume
    cbf_ts
       nd array of 3D or 4D computed cbf
    gm,wm,csf 
       numpy array of grey matter, whitematter, and csf
    mask 
       numpy array of mask 

    reference:

    """
    gm[gm < thresh] = 0
    gm[gm > 0] = 1
    wm[wm < thresh] = 0
    wm[wm > 0] = 1
    csf[csf < thresh] = 0
    csf[csf > 0] = 1
    # get the total number of voxle within csf,gm and wm
    nogm = np.sum(gm == 1) - 1
    nowm = np.sum(wm == 1) - 1
    nocf = np.sum(csf == 1) - 1
    mask1 = gm + wm + csf
    # msk=sum(mask>0)
    # mean  of times series cbf within greymatter
    mgmts = np.squeeze(np.mean(cbfts[gm == 1, :], axis=0))
    # robiust mean and meadian
    from scipy.stats import median_absolute_deviation
    medmngm = np.median(mgmts)
    sdmngm = median_absolute_deviation(mgmts)/0.675
    indx = 1*(np.abs(mgmts-medmngm) > (2.5*sdmngm))
    R = np.mean(cbfts[:, :, :, indx == 0], axis=3)
    V = nogm*np.var(R[gm == 1]) + nowm*np.var(R[wm == 1]) + nocf*np.var(R[csf == 1])
    V1 = V+1
    while V < V1:
        V1 = V
        CC = np.zeros(cbfts.shape[3])*(-2)
        for s in range(cbfts.shape[3]):
            if indx[s] != 0:
                break
            else:
                tmp1 = cbfts[:, :, :, s]
                CC[s] = np.corrcoef(R[mask1 > 0], tmp1[mask1 > 0])[0][1]
        inx = np.argmax(CC)
        indx[inx] = 2
        R = np.mean(cbfts[:, :, :, indx == 0], axis=3)
        V = nogm*np.var(R[gm == 1]) + nowm*np.var(R[wm == 1]) + nocf*np.var(R[csf == 1])
    cbfts_recon = cbfts[:, :, :, indx == 0]
    cbfts_recon1 = np.zeros_like(cbfts_recon)
    for i in range(cbfts_recon.shape[3]):
        cbfts_recon1[:, :, :, i] = cbfts_recon[:, :, :, i]*mask
    cbfts_recon1 = np.nan_to_num(cbfts_recon1)
    return cbfts_recon1, indx


def _roubustfit(Y, mu, Globalprior, modrobprior, lmd=0, localprior=0, wfun='huber', tune=1.345,
                flagstd=1, flagmodrobust=1, flagprior=1, thresh=0.7):
    """
    robust fit 
    """
    dimcbf = Y.shape
    priow = np.ones([dimcbf[0], dimcbf[1]])
    sw = 1
    X = priow
    b = (np.sum(X*Y, axis=0)+mu*Globalprior+lmd*localprior)/(np.sum(X*X, axis=0)+mu+lmd)
    b0 = np.repeat(0, len(b))
    h1 = X/np.power(np.tile(np.sqrt(np.sum(X*X, axis=0)), (dimcbf[0], 1)), 2)
    h0 = 0.9999*np.ones([dimcbf[0], dimcbf[1]])
    h = np.minimum(h0, h1)
    adjfactor = 1/(np.sqrt(1-h/priow))
    tiny_s = (1e-6)*(np.std(h, axis=0))
    tiny_s[tiny_s == 0] = 1
    D = np.sqrt(np.finfo(float).eps)
    iter = 0
    interlim = 10
    while iter < interlim:
        print('iteration  ', iter, "\n")
        iter = iter + 1
        check1 = np.subtract(np.abs(b-b0), (D*np.maximum(np.abs(b), np.abs(b0))))
        check1[check1 > 0] = 0
        if any(check1):
            print(' \n converged after ', iter, "iterations\n")
            break
        r = Y - X*(np.tile(b, (dimcbf[0], 1)))
        radj = r * adjfactor/sw
        if flagstd == 1:
            s = np.sqrt(np.mean(np.power(radj, 2), axis=0))
        else:
            rs = np.sort(np.abs(radj), axis=0)
            s = np.median(rs, axis=0)/0.6745
        rx1 = radj*(1-flagmodrobust*np.exp(-np.tile(modrobprior, (dimcbf[0], 1))))
        rx2 = np.tile(np.maximum(s, tiny_s)*tune, (dimcbf[0], 1))
        r1 = rx1 / rx2
        w, _ = _weightfun(r1, wfun)
        b0 = b
        z = np.sqrt(w)
        x = X*z
        yz = Y*z
        b = (np.sum(x*yz, axis=0)+mu*Globalprior+lmd*localprior)/(np.sum(x*x, axis=0)+mu+lmd)
        b = np.nan_to_num(b)
    return b


def _scrubcbf(cbf_ts, gm, wm, csf, mask, wfun='huber', thresh=0.7):
    
    """ 
    scrub algorithms by Sudipto
    cbf_ts
       nd array of 3D or 4D computed cbf
    gm,wm,csf 
       numpy array of grey matter, whitematter, and csf
    mask 
       numpy array of mask 
    
    wf 
      wave function

    reference:

    """

    gm = mask*gm
    wm = mask*wm
    csf = csf*mask
    gmidx = gm[mask == 1]
    gmidx[gmidx < thresh] = 0
    gmidx[gmidx > 0] = 1
    wmidx = wm[mask == 1]
    wmidx[wmidx < thresh] = 0
    wmidx[wmidx > 0] = 1
    csfidx = csf[mask == 1]
    csfidx[csfidx < thresh] = 0
    csfidx[csfidx > 0] = 1
    # midx = mask[mask==1]
    meancbf = np.mean(cbf_ts, axis=3)
    y = np.transpose(cbf_ts[mask == 1, :, ])
    VV = np.var(y, axis=0)
    thresh1, thresh3 = _getchisquare(y.shape[0])
    mu1 = VV/(np.median(VV[gmidx == 1])*thresh3)
    mu = ((mu1 > thresh1) & (mu1 < 10*thresh1))*(mu1-thresh1) \
        + (mu1 >= 10*thresh1)*(1/(2*thresh1*10) * np.power(mu1, 2))+(thresh1*10/2 - thresh1)
    M = meancbf*mask
    M[mask == 1] = mu
    modrobprior = mu/10
    gmidx2 = 1*([gm.flatten() > thresh] and [M.flatten() == 0]
                and [wm.flatten() > csf.flatten()])[0]
    wmidx2 = 1*([wm.flatten() > thresh] and [M.flatten() == 0]
                and [gm.flatten() > csf.flatten()])[0]
    if np.sum(gmidx2) == 0 or np.sum(wmidx2) == 0:
        gmidx2 = 1*(gm.flatten() > thresh)
        wmidx2 = 1*(wm.flatten() > thresh)
    idxx = gmidx2 + wmidx2
    idxx[idxx > 0] = 1
    X = np.zeros([len(idxx), 2])
    X[:, 0] = gm.flatten()[gm.flatten() >= (0)]*idxx
    X[:, 1] = wm.flatten()[wm.flatten() >= (0)]*idxx
    A = (meancbf.flatten()[idxx >= 0])*idxx
    c = np.linalg.lstsq(X, A)[0]
    Globalpriorfull = c[0]*gm.flatten() + c[1]*wm.flatten()
    Globalprior = Globalpriorfull[mask.flatten() == 1]
    localprior = 0
    lmd = 0
    tune = _tune(wfun=wfun)
    bb = _roubustfit(Y=y, mu=mu, Globalprior=Globalprior, modrobprior=modrobprior,
                     lmd=lmd, localprior=localprior, wfun=wfun, tune=tune, flagstd=1,
                     flagmodrobust=1, flagprior=1, thresh=0.7)
    newcbf = meancbf*mask
    newcbf[mask == 1] = bb
    newcbf=np.nan_to_num(newcbf)
    return newcbf

# basil and pvcorr


class _BASILCBFInputSpec(FSLCommandInputSpec):
    # We use position args here as list indices - so a negative number
    # will put something on the end
    in_file = File(
        exists=True,
        desc="input file cbf after substracting tag-control or control-tag",
        argstr=" -i %s",
        position=0,
        mandatory=True,
    )
    mask = File(exists=True, argstr=" -m %s ", desc="mask in the same space as in_infile",
                mandatory=True,)
    mzero = File(exists=True, argstr=" -c %s ", desc='m0 scan', mandatory=False)
    m0scale = traits.Float(desc='calibration of asl', argstr=" --cgain %.2f ", mandatory=True)
    m0tr = traits.Float(desc='Mzero TR', argstr=" --tr %.2f ", mandatory=True,)
    tis = traits.Str(desc='recovery time =plds+bolus', argstr=" --tis %s ", mandatory=True,)
    pcasl = traits.Bool(desc='label type:defualt is PASL', argstr=" --casl ",
                        mandatory=False, default_value=False)
    bolus = traits.Float(desc='bolus or tau: label duration', argstr=" --bolus %.2f ",
                         mandatory=True)
    pvc = traits.Bool(desc='calibration of asl', mandatory=False, argstr=" --pvcorr ",
                      default_value=True)
    pvgm = File(exists=True, mandatory=False, desc='grey matter probablity matter ',
                argstr=" --pvgm %s ")
    pvwm = File(exists=True, mandatory=False, desc='white matter probablity matter ',
                argstr=" --pvwm %s ")
    out_basename = File(desc="base name of output files", argstr=" -o %s ", mandatory=True)
    out_cbfb = File(exists=False, desc='cbf with spatial correction')
    out_cbfpv = File(exists=False, desc='cbf with spatial correction')
    out_att = File(exists=False, desc='aretrial transist time')
    # environ=traits.Str('FSLOUTPUTTYPE': 'NIFTI_GZ'}


class _BASILCBFOutputSpec(TraitedSpec):
    out_cbfb = File(exists=False, desc='cbf with spatial correction')
    out_cbfpv = File(exists=False, desc='cbf with spatial correction')
    out_att = File(exists=False, desc='aretrial transist time')


class BASILCBF(FSLCommand):
    r"""
    oxford asl 
    https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/BASIL
    """
    _cmd = " oxford_asl "
    input_spec = _BASILCBFInputSpec
    output_spec = _BASILCBFOutputSpec

    def _run_interface(self, runtime):
        runtime = super(BASILCBF, self)._run_interface(runtime)
        return runtime

    def _gen_outfilename(self, suffix):
        if isdefined(self.inputs.in_file):
            out_file = self._gen_fname(self.inputs.in_file, suffix=suffix)
        return os.path.abspath(out_file)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        # outputs["out_cbfb"]=self.inputs.out_basename+'/basilcbf.nii.gz'
        outputs["out_cbfb"] = fname_presuffix(self.inputs.mask, suffix='_cbfbasil')
        from shutil import copyfile
        copyfile(self.inputs.out_basename+'/native_space/perfusion_calib.nii.gz',
                 outputs["out_cbfb"])
        
        
            # outputs["out_att"]=self.inputs.out_basename+'/arrivaltime.nii.gz'
        outputs["out_att"] = fname_presuffix(self.inputs.mask, suffix='_arrivaltime')
        copyfile(self.inputs.out_basename+'/native_space/arrival.nii.gz', outputs["out_att"])
        self.inputs.out_att = os.path.abspath(outputs["out_att"])
        
        # outputs["out_cbfpv"]=self.inputs.out_basename+'/basilcbfpv.nii.gz'
        outputs["out_cbfpv"] = fname_presuffix(self.inputs.mask, suffix='_cbfbasilpv')
        copyfile(self.inputs.out_basename+'/native_space/pvcorr/perfusion_calib.nii.gz',
                 outputs["out_cbfpv"])
        self.inputs.out_cbfb = os.path.abspath(outputs["out_cbfb"])
        self.inputs.out_cbfpv = os.path.abspath(outputs["out_cbfpv"])
        return outputs


class _qccbfInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='original asl_file')
    in_meancbf = File(exists=True, mandatory=True, desc='cbf img')
    in_avgscore = File(exists=True, mandatory=False, desc='cbf img')
    in_scrub = File(exists=True, mandatory=False, desc='cbf img')
    in_basil = File(exists=True, mandatory=False, desc='cbf img')
    in_pvc = File(exists=True, mandatory=False, desc='cbf img')
    in_greyM = File(exists=True, mandatory=True, desc='grey  matter')
    in_whiteM = File(exists=True, mandatory=True, desc='white  matter')
    in_csf = File(exists=True, mandatory=True, desc='csf')
    in_confmat = File(exists=True, mandatory=False, desc=' cofnound matrix')
    in_aslmask = File(exists=True, mandatory=True, desc='asl mask in native space')
    in_t1mask = File(exists=True, mandatory=True, desc='t1wmask in native space ')
    in_aslmaskstd = File(exists=True, mandatory=False, desc='asl mask in native space')
    in_templatemask = File(exists=True, mandatory=False, desc='template mask or image')
    qc_file = File(exists=False, mandatory=False, desc='qc file ')


class _qccbfOutputSpec(TraitedSpec):
    qc_file = File(exists=False, desc='qc file ')


class qccbf(SimpleInterface):
    r""""
     compute qc from confound regressors 
     and cbf maps, 
     coregistration and regsitration indexes

    """

    input_spec = _qccbfInputSpec
    output_spec = _qccbfOutputSpec

    def _run_interface(self, runtime):

        time1 = pd.read_csv(self.inputs.in_confmat, sep='\t')
        time1.fillna(0, inplace=True)
        fd = np.mean(time1['framewise_displacement'])
        rms = time1[['rot_x', 'rot_y', 'rot_z']]
        rms1 = rms.pow(2)
        rms = np.mean(np.sqrt(rms1.sum(axis=1)/3))
        regDC = dc(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regJC = jc(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regCC = crosscorr(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regCov = coverage(self.inputs.in_aslmask, self.inputs.in_t1mask)

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            normDC = dc(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normJC = jc(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normCC = crosscorr(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normCov = coverage(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)

        meancbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                              csf=self.inputs.in_csf, img=self.inputs.in_meancbf, thresh=0.7)
        meancbf = globalcbf(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                            csf=self.inputs.in_csf, cbf=self.inputs.in_meancbf, thresh=0.7)

        if self.inputs.in_avgscore:
            scorecbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                               csf=self.inputs.in_csf, img=self.inputs.in_avgscore, thresh=0.7)
            scrub_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                            csf=self.inputs.in_csf, img=self.inputs.in_scrub, thresh=0.7)
            negscore = negativevoxel(cbf=self.inputs.in_avgscore, gm=self.inputs.in_greyM, thresh=0.7)
            negscrub = negativevoxel(cbf=self.inputs.in_scrub, gm=self.inputs.in_greyM, thresh=0.7)
        else:
            scorecbf_qei = 0
            scrub_qei = 0 
            negscore = 0
            negscrub = 0

        if self.inputs.in_basil:
            basilcbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                               csf=self.inputs.in_csf, img=self.inputs.in_basil, thresh=0.7)
            pvcbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                            csf=self.inputs.in_csf, img=self.inputs.in_pvc, thresh=0.7)
            negbasil = negativevoxel(cbf=self.inputs.in_basil, gm=self.inputs.in_greyM, thresh=0.7)
            negpvc = negativevoxel(cbf=self.inputs.in_pvc, gm=self.inputs.in_greyM, thresh=0.7)
        else:
            basilcbf_qei = 0
            pvcbf_qei = 0 
            negbasil = 0
            negpvc = 0
        
        
        gwratio = np.divide(meancbf[0], meancbf[1])
        negcbf = negativevoxel(cbf=self.inputs.in_meancbf, gm=self.inputs.in_greyM, thresh=0.7)
        
        

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            dict1 = {'FD': [fd], 'relRMS': [rms], 'coregDC': [regDC], 'coregJC': [regJC],
                     'coregCC': [regCC], 'coregCOV': [regCov], 'normDC': [normDC],
                     'normJC': [normJC], 'normCC': [normCC], 'normCOV': [normCov],
                     'cbfQEI': [meancbf_qei], 'scoreQEI': [scorecbf_qei], 'scrubQEI': [scrub_qei],
                     'basilQEI': [basilcbf_qei], 'pvcQEI': [pvcbf_qei], 'GMmeanCBF': [meancbf[0]],
                     'WMmeanCBF': [meancbf[1]], 'Gm_Wm_CBF_ratio': [gwratio],
                     'NEG_CBF_PERC': [negcbf], 'NEG_SCORE_PERC': [negscore],
                     'NEG_SCRUB_PERC': [negscrub], 'NEG_BASIL_PERC': [negbasil],
                     'NEG_PVC_PERC': [negpvc]}
        else:
            dict1 = {'FD': [fd], 'relRMS': [rms], 'coregDC': [regDC], 'coregJC': [regJC],
                     'coregCC': [regCC], 'coregCOV': [regCov],
                     'cbfQEI': [meancbf_qei], 'scoreQEI': [scorecbf_qei], 'scrubQEI': [scrub_qei],
                     'basilQEI': [basilcbf_qei], 'pvcQEI': [pvcbf_qei], 'GMmeanCBF': [meancbf[0]],
                     'WMmeanCBF': [meancbf[1]], 'Gm_Wm_CBF_ratio': [gwratio],
                     'NEG_CBF_PERC': [negcbf], 'NEG_SCORE_PERC': [negscore],
                     'NEG_SCRUB_PERC': [negscrub], 'NEG_BASIL_PERC': [negbasil],
                     'NEG_PVC_PERC': [negpvc]}
        _, file1 = os.path.split(self.inputs.in_file)
        bb = file1.split('_')
        dict2 = {}
        for i in range(len(bb)-1):
            dict2.update({bb[i].split('-')[0]: bb[i].split('-')[1]})
        dict2.update(dict1)

        df = pd.DataFrame(dict2)

        self._results['qc_file'] = fname_presuffix(self.inputs.in_meancbf, suffix='qc_cbf.csv',
                                                   newpath=runtime.cwd, use_ext=False)
        df.to_csv(self._results['qc_file'], index=False, header=True)

        self.inputs.qc_file = os.path.abspath(self._results['qc_file'])
        return runtime
        

class _qccbfgeInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='original asl_file')
    in_meancbf = File(exists=True, mandatory=True, desc='cbf img')
    in_avgscore = File(exists=True, mandatory=False, desc='cbf img')
    in_scrub = File(exists=True, mandatory=False, desc='cbf img')
    in_basil = File(exists=True, mandatory=False, desc='cbf img')
    in_pvc = File(exists=True, mandatory=False, desc='cbf img')
    in_greyM = File(exists=True, mandatory=True, desc='grey  matter')
    in_whiteM = File(exists=True, mandatory=True, desc='white  matter')
    in_csf = File(exists=True, mandatory=True, desc='csf')
    in_aslmask = File(exists=True, mandatory=True, desc='asl mask in native space')
    in_t1mask = File(exists=True, mandatory=True, desc='t1wmask in native space ')
    in_aslmaskstd = File(exists=True, mandatory=False, desc='asl mask in native space')
    in_templatemask = File(exists=True, mandatory=False, desc='template mask or image')
    qc_file = File(exists=False, mandatory=False, desc='qc file ')


class _qccbfgeOutputSpec(TraitedSpec):
    qc_file = File(exists=False, desc='qc file ')


class qccbfge(SimpleInterface):
    r""""
     compute qc from confound regressors 
     and cbf maps, 
     coregistration and regsitration indexes

    """

    input_spec = _qccbfInputSpec
    output_spec = _qccbfOutputSpec

    def _run_interface(self, runtime):
        regDC = dc(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regJC = jc(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regCC = crosscorr(self.inputs.in_aslmask, self.inputs.in_t1mask)
        regCov = coverage(self.inputs.in_aslmask, self.inputs.in_t1mask)

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            normDC = dc(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normJC = jc(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normCC = crosscorr(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)
            normCov = coverage(self.inputs.in_aslmaskstd, self.inputs.in_templatemask)

        meancbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                              csf=self.inputs.in_csf, img=self.inputs.in_meancbf, thresh=0.7)
        meancbf = globalcbf(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                            csf=self.inputs.in_csf, cbf=self.inputs.in_meancbf, thresh=0.7)
        

        if self.inputs.in_avgscore:
            scorecbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                               csf=self.inputs.in_csf, img=self.inputs.in_avgscore, thresh=0.7)
            scrub_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                            csf=self.inputs.in_csf, img=self.inputs.in_scrub, thresh=0.7)
            negscore = negativevoxel(cbf=self.inputs.in_avgscore, gm=self.inputs.in_greyM, thresh=0.7)
            negscrub = negativevoxel(cbf=self.inputs.in_scrub, gm=self.inputs.in_greyM, thresh=0.7)
        else:
            scorecbf_qei = 0
            scrub_qei = 0 
            negscore = 0
            negscrub = 0

        if self.inputs.in_basil:
            basilcbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                               csf=self.inputs.in_csf, img=self.inputs.in_basil, thresh=0.7)
            pvcbf_qei = cbf_qei(gm=self.inputs.in_greyM, wm=self.inputs.in_whiteM,
                            csf=self.inputs.in_csf, img=self.inputs.in_pvc, thresh=0.7)
            negbasil = negativevoxel(cbf=self.inputs.in_basil, gm=self.inputs.in_greyM, thresh=0.7)
            negpvc = negativevoxel(cbf=self.inputs.in_pvc, gm=self.inputs.in_greyM, thresh=0.7)
        else:
            basilcbf_qei = 0
            pvcbf_qei = 0 
            negbasil = 0
            negpvc = 0
        gwratio = np.divide(meancbf[0], meancbf[1])
        negcbf = negativevoxel(cbf=self.inputs.in_meancbf, gm=self.inputs.in_greyM, thresh=0.7)

        if self.inputs.in_aslmaskstd and self.inputs.in_templatemask:
            dict1 = {'FD': 0, 'relRMS': 0, 'coregDC': [regDC], 'coregJC': [regJC],
                     'coregCC': [regCC], 'coregCOV': [regCov], 'normDC': [normDC],
                     'normJC': [normJC], 'normCC': [normCC], 'normCOV': [normCov],
                     'cbfQEI': [meancbf_qei], 'scoreQEI': [scorecbf_qei], 'scrubQEI': [scrub_qei],
                     'basilQEI': [basilcbf_qei], 'pvcQEI': [pvcbf_qei], 'GMmeanCBF': [meancbf[0]],
                     'WMmeanCBF': [meancbf[1]], 'Gm_Wm_CBF_ratio': [gwratio],
                     'NEG_CBF_PERC': [negcbf], 'NEG_SCORE_PERC': [negscore],
                     'NEG_SCRUB_PERC': [negscrub], 'NEG_BASIL_PERC': [negbasil],
                     'NEG_PVC_PERC': [negpvc]}
        else:
            dict1 = {'FD': 0, 'relRMS': 0, 'coregDC': [regDC], 'coregJC': [regJC],
                     'coregCC': [regCC], 'coregCOV': [regCov],
                     'cbfQEI': [meancbf_qei], 'scoreQEI': [scorecbf_qei], 'scrubQEI': [scrub_qei],
                     'basilQEI': [basilcbf_qei], 'pvcQEI': [pvcbf_qei], 'GMmeanCBF': [meancbf[0]],
                     'WMmeanCBF': [meancbf[1]], 'Gm_Wm_CBF_ratio': [gwratio],
                     'NEG_CBF_PERC': [negcbf], 'NEG_SCORE_PERC': [negscore],
                     'NEG_SCRUB_PERC': [negscrub], 'NEG_BASIL_PERC': [negbasil],
                     'NEG_PVC_PERC': [negpvc]}
        _, file1 = os.path.split(self.inputs.in_file)
        bb = file1.split('_')
        dict2 = {}
        for i in range(len(bb)-1):
            dict2.update({bb[i].split('-')[0]: bb[i].split('-')[1]})
        dict2.update(dict1)

        df = pd.DataFrame(dict2)

        self._results['qc_file'] = fname_presuffix(self.inputs.in_meancbf, suffix='qc_cbf.csv',
                                                   newpath=runtime.cwd, use_ext=False)
        df.to_csv(self._results['qc_file'], index=False, header=True)

        self.inputs.qc_file = os.path.abspath(self._results['qc_file'])
        return runtime

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
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
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
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))

    intersection = np.count_nonzero(input1 & input2)
    union = np.count_nonzero(input1 | input2)

    jc = float(intersection) / float(union)

    return jc


def crosscorr(input1, input2):
    r"""
    cross correlation
    computer compute cross correction bewteen input mask
    """
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool)).flatten()
    input2 = np.atleast_1d(input2.astype(np.bool)).flatten()
    cc = np.corrcoef(input1, input2)[0][1]
    return cc


def coverage(input1, input2):
    """
    estimate the coverage between  two mask
    """
    input1 = nb.load(input1).get_fdata()
    input2 = nb.load(input2).get_fdata()
    input1 = np.atleast_1d(input1.astype(np.bool))
    input2 = np.atleast_1d(input2.astype(np.bool))
    intsec = np.count_nonzero(input1 & input2)
    if np.sum(input1) > np.sum(input2):
        smallv = np.sum(input2)
    else:
        smallv = np.sum(input1)
    cov = float(intsec)/float(smallv)
    return cov


def globalcbf(cbf, gm, wm, csf, thresh=0.7):
    cbf = nb.load(cbf).get_fdata()
    gm = nb.load(gm).get_fdata()
    wm = nb.load(wm).get_fdata()
    csf = nb.load(csf).get_fdata()
    b1 = [gm < thresh]
    gm[b1] = 0
    bx = cbf[gm > 0]
    b2 = [wm < thresh]
    wm[b2] = 0
    by = cbf[wm > 0]
    b3 = [csf < thresh]
    csf[b3] = 0
    bz = cbf[csf > 0]
    return np.mean(bx), np.mean(by), np.mean(bz)


def cbf_qei(gm, wm, csf, img, thresh=0.7):
    """
    Quality evaluation index of CBF base on Sudipto Dolui work 
    Dolui S., Wolf R. & Nabavizadeh S., David W., Detre, J. (2017). 
    Automated Quality Evaluation Index for 2D ASL CBF Maps. ISMR 2017

    """
    def fun1(x, xdata):
        d1 = np.exp(-(x[0])*np.power(xdata, x[1]))
        return(d1)

    def fun2(x, xdata):
        d1 = 1-np.exp(-(x[0])*np.power(xdata, x[1]))
        return(d1)

    x1 = [0.054, 0.9272]
    x2 = [2.8478, 0.5196]
    x4 = [3.0126, 2.4419]
    scbf = smooth_image(nb.load(img), fwhm=5).get_fdata()
    if len(scbf.shape) > 3:
        scbf = scbf[:, :, :, 0]
    # load prob maps
    gmm = nb.load(gm).get_fdata()
    wmm = nb.load(wm).get_fdata()
    ccf = nb.load(csf).get_fdata()
    if len(gmm.shape) > 3:
        gmm = gmm[:, :, :, 0]
        wmm = wmm[:, :, :, 0]
        ccf = ccf[:, :, :, 0]
    pbcf = 2.5*gmm+wmm  # gmm is 2.5 times wm
    msk = np.array((scbf != 0) & (scbf != np.nan) & (pbcf != np.nan)).astype(int)

    gm1 = np.array(gmm > thresh)
    wm1 = np.array(wmm > thresh)
    cc1 = np.array(ccf > thresh)
    r1 = np.array([0, np.corrcoef(scbf[msk == 1], pbcf[msk == 1])[1, 0]]).max()

    V = ((np.sum(gm1)-1)*np.var(scbf[gm1 > 0])+(np.sum(wm1)-1)*np.var(scbf[wm1 > 0])
         + (np.sum(cc1)-1) * np.var(scbf[cc1 > 0]))/(np.sum(gm1 > 0)+np.sum(wm1 > 0)
                                                     + np.sum(cc1 > 0)-3)

    negGM = np.sum(scbf[gm1] < 0)/(np.sum(gm1))
    GMCBF = np.mean(scbf[gm1])
    CV = V/np.abs(GMCBF)
    Q = [fun1(x1, CV), fun1(x2, negGM), fun2(x4, r1)]
    return gmean(Q)


def negativevoxel(cbf, gm, thresh=0.7):
    """
    percentage of negative voxel within
    grey matter voxel
    """
    gm = nb.load(gm).get_fdata()
    cbf = nb.load(cbf).get_fdata()
    gm1 = np.array(gm > thresh)
    gm1[gm1 > 0] = 1
    npgm = np.sum(gm1)
    cbfgm = np.array(cbf < 0)
    cbfgm[cbfgm < 0] = 1
    ncbfgm = np.sum(np.multiply(cbfgm, gm1))
    pernegcbf = np.multiply(np.divide(ncbfgm, npgm), 100)
    return pernegcbf


def get_atlas(atlasname):
    if atlasname == 'HarvardOxford':
        atlasfile = pkgrf('aslprep', 'data/atlas/HarvardOxford/HarvardOxfordMNI.nii.gz')
        atlasdata = pkgrf('aslprep', 'data/atlas/HarvardOxford/HarvardOxfordNodeNames.txt')
        atlaslabel = pkgrf('aslprep', 'data/atlas/HarvardOxford/HarvardOxfordNodeIndex.1D')
    elif atlasname == 'schaefer200x7':
        atlasfile = pkgrf('aslprep', 'data/atlas/schaefer200x7/schaefer200x7MNI.nii.gz')
        atlasdata = pkgrf('aslprep', 'data/atlas/schaefer200x7/schaefer200x7NodeNames.txt')
        atlaslabel = pkgrf('aslprep', 'data/atlas/schaefer200x7/schaefer200x7NodeIndex.1D')
    elif atlasname == 'schaefer200x17':
        atlasfile = pkgrf('aslprep', 'data/atlas/schaefer200x17/schaefer200x17MNI.nii.gz')
        atlasdata = pkgrf('aslprep', 'data/atlas/schaefer200x17/schaefer200x17NodeNames.txt')
        atlaslabel = pkgrf('aslprep', 'data/atlas/schaefer200x17/schaefer200x17NodeIndex.1D')
    elif atlasname == 'schaefer400x7':
        atlasfile = pkgrf('aslprep', 'data/atlas/schaefer400x7/schaefer400x7MNI.nii.gz')
        atlasdata = pkgrf('aslprep', 'data/atlas/schaefer400x7/schaefer400x7NodeNames.txt')
        atlaslabel = pkgrf('aslprep', 'data/atlas/schaefer200x17/schaefer200x17NodeIndex.1D')
    elif atlasname == 'schaefer400x17':
        atlasfile = pkgrf('aslprep', 'data/atlas/schaefer400x17/schaefer400x17MNI.nii.gz')
        atlasdata = pkgrf('aslprep', 'data/atlas/schaefer400x17/schaefer400x17NodeNames.txt')
        atlaslabel = pkgrf('aslprep', 'data/atlas/schaefer400x17/schaefer400x17NodeIndex.1D')
    else:
        raise RuntimeError('atlas not available')
    return atlasfile, atlasdata, atlaslabel


def cbfroiquant(roi_file, roi_label, cbfmap):
    data = nb.load(cbfmap).get_data()
    roi = nb.load(roi_file).get_data()
    roi_labels = np.loadtxt(roi_label)
    if (data.shape != roi.shape):
        raise ValueError("Image-shapes do not match")
    # if roi_labels is None:
    # roi_labels = np.unique(roi)
    mean_vals = []
    for roi_label in roi_labels:
        mean_vals.append(np.mean(data[roi == roi_label]))
    return mean_vals


class _cbfroiquantInputSpec(BaseInterfaceInputSpec):
    in_cbf = File(exists=True, mandatory=True, desc='cbf img')
    atlasfile = File(exists=True, mandatory=True, desc='data')
    atlasdata = File(exists=True, mandatory=True, desc='data')
    atlaslabel = File(exists=True, mandatory=True, desc='data')
    atlascsv = File(exists=False, mandatory=False, desc='harvard output csv')


class _cbfroiquantOutputSpec(TraitedSpec):
    atlascsv = File(exists=False, desc='harvard output csv')


class cbfqroiquant(SimpleInterface):
    input_spec = _cbfroiquantInputSpec
    output_spec = _cbfroiquantOutputSpec

    def _run_interface(self, runtime):

        self._results['atlascsv'] = fname_presuffix(self.inputs.in_cbf, suffix='atlas.csv',
                                                    newpath=runtime.cwd, use_ext=False)
        roiquant = cbfroiquant(roi_label=self.inputs.atlaslabel, roi_file=self.inputs.atlasfile,
                               cbfmap=self.inputs.in_cbf)
        data1 = pd.read_table(self.inputs.atlasdata, header=None, index_col=None, sep='\t')
        bb = list(data1.values.tolist())
        flattened = [val for sublist in bb for val in sublist]
        datat = pd.DataFrame([flattened, roiquant])
        datat.to_csv(self._results['atlascsv'], header=None, index=None)
        return runtime


def regmotoasl(asl,m0file,m02asl):
    from nipype.interfaces import fsl
    meanasl = fsl.MeanImage(); meanasl.inputs.in_file = asl
    meanasl.inputs.out_file = fname_presuffix(asl,suffix='_meanasl')
    meanasl.run()
    meanm0 = fsl.MeanImage(); meanm0.inputs.in_file = m0file
    meanm0.inputs.out_file = fname_presuffix(asl,suffix='_meanm0')
    meanm0.run()
    flt = fsl.FLIRT(bins=640, cost_func='mutualinfo')
    flt.inputs.in_file = meanm0.inputs.out_file 
    flt.inputs.reference = meanasl.inputs.out_file
    flt.inputs.out_file = m02asl
    flt.run()
    return m02asl    


def readjson(jsonfile):
    import json
    with open(jsonfile) as f:
        data = json.load(f)
    return data

class _extractCBInputSpec(BaseInterfaceInputSpec):
    in_asl = File(exists=True, mandatory=True, desc='raw asl file')
    in_aslmask = File(exists=True, mandatory=True, desct='asl mask')
    file_type = traits.Str(desc='file type, c for cbf, d for deltam',mandatory=True)
    out_file = File(exists=False, mandatory=False, desc='cbf or deltam')


class _extractCBOutputSpec(TraitedSpec):
    out_file= File(exists=False, desc='cbf or deltalm')
    


class extractCB(SimpleInterface):
    r"""
    the code refine the asl mask with t1w mask
    the output is refined asl mask

    """
    input_spec = _extractCBInputSpec
    output_spec = _extractCBOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = fname_presuffix(self.inputs.in_aslmask,
                                                   suffix='_cbfdeltam', newpath=runtime.cwd)
        filex=self.inputs.in_asl
        aslcontext = pd.read_csv(filex.replace('_asl.nii.gz', '_aslcontext.tsv'))
        idasl = aslcontext['volume_type'].tolist()
        
        fdata=nb.load(filex).get_fdata()
        img=nb.load(filex)
        if self.inputs.file_type == 'd':
            dlist = [i for i in range(0, len(idasl)) if idasl[i] == 'deltam']
        elif self.inputs.file_type == 'c':
            dlist = [i for i in range(0, len(idasl)) if idasl[i] == 'CBF']
        
        if len(fdata.shape) < 4:
            newdata = nb.Nifti1Image(dataobj=fdata,affine=img.affine,header=img.header)
        else:
            ffdata=fdata[:, :, :, dlist]
            newdata = nb.Nifti1Image(dataobj=ffdata,affine=img.affine,header=img.header)
        
        newdata.to_filename(self._results['out_file'])

        return runtime