import nibabel as nb
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import gridspec as mgs
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.colorbar import ColorbarBase

from nilearn.plotting import plot_img
from nilearn.signal import clean
from nilearn.image import resample_to_img
from nilearn._utils import check_niimg_4d
from nilearn._utils.niimg import _safe_get_data


def plot_carpet(img, atlaslabels, detrend=True, nskip=0, size=(950, 800),
                subplot=None, title=None, legend=False):
    """
    Plot an image representation of voxel intensities across time also know
    as the "carpet plot" or "Power plot". See Jonathan Power Neuroimage
    2017 Jul 1; 154:150-158.
    Parameters
    ----------
        img : Niimg-like object
            See http://nilearn.github.io/manipulating_images/input_output.html
            4D input image
        atlaslabels: ndarray
            A 3D array of integer labels from an atlas, resampled into ``img`` space.
        detrend : boolean, optional
            Detrend and standardize the data prior to plotting.
        nskip : int
            Number of volumes at the beginning of the scan marked as nonsteady state.
        axes : matplotlib axes, optional
            The axes used to display the plot. If None, the complete
            figure is used.
        title : string, optional
            The title displayed on the figure.
        output_file : string, or None, optional
            The name of an image file to export the plot to. Valid extensions
            are .png, .pdf, .svg. If output_file is not None, the plot
            is saved to a file, and the display is closed.
        legend : bool
            Whether to render the average functional series with ``atlaslabels`` as
            overlay.
    """
    img_nii = check_niimg_4d(img, dtype='auto')
    func_data = _safe_get_data(img_nii, ensure_finite=True)

    # Define TR and number of frames
    tr = img_nii.header.get_zooms()[-1]
    ntsteps = func_data.shape[-1]

    data = func_data[atlaslabels > 0].reshape(-1, ntsteps)
    seg = atlaslabels[atlaslabels > 0].reshape(-1)

    # Apply lookup table
    lut = np.unique(atlaslabels)
    newsegm = lut[seg.astype(int)]

    p_dec = 1 + data.shape[0] // size[0]
    if p_dec:
        data = data[::p_dec, :]
        newsegm = newsegm[::p_dec]

    t_dec = 1 + data.shape[1] // size[1]
    if t_dec:
        data = data[:, ::t_dec]

    # Detrend data
    v = (None, None)
    if detrend:
        data = clean(data.T, t_r=tr).T
        v = (-2, 2)

    # Order following segmentation labels
    order = np.argsort(newsegm)[::-1]

    # If subplot is not defined
    if subplot is None:
        subplot = mgs.GridSpec(1, 1)[0]

    # Define nested GridSpec
    wratios = [1, 100, 20]
    gs = mgs.GridSpecFromSubplotSpec(1, 2 + int(legend), subplot_spec=subplot,
                                     width_ratios=wratios[:2 + int(legend)],
                                     wspace=0.0)

    mycolors = ListedColormap(cm.get_cmap('tab10').colors[:4][::-1])

    # Segmentation colorbar
    ax0 = plt.subplot(gs[0])
    ax0.set_yticks([])
    ax0.set_xticks([])
    ax0.imshow(newsegm[order, np.newaxis], interpolation='none', aspect='auto',
               cmap=mycolors, vmin=1, vmax=4)
    ax0.grid(False)
    ax0.spines["left"].set_visible(False)
    ax0.spines["bottom"].set_color('none')
    ax0.spines["bottom"].set_visible(False)

    print(data.shape)
    print(data[order].shape)

    # Carpet plot
    ax1 = plt.subplot(gs[1])
    ax1.imshow(data[order, ...], interpolation='nearest', aspect='auto', cmap='gray',
               vmin=v[0], vmax=v[1])

    ax1.grid(False)
    ax1.set_yticks([])
    ax1.set_yticklabels([])

    # Set 10 frame markers in X axis
    interval = max((int(data.shape[-1] + 1) // 10, int(data.shape[-1] + 1) // 5, 1))
    xticks = list(range(0, data.shape[-1])[::interval])
    ax1.set_xticks(xticks)
    ax1.set_xlabel('time (s)')
    labels = tr * (np.array(xticks)) * t_dec
    ax1.set_xticklabels(['%.02f' % t for t in labels.tolist()], fontsize=5)

    # Remove and redefine spines
    for side in ["top", "right"]:
        # Toggle the spine objects
        ax0.spines[side].set_color('none')
        ax0.spines[side].set_visible(False)
        ax1.spines[side].set_color('none')
        ax1.spines[side].set_visible(False)

    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.spines["bottom"].set_visible(False)
    ax1.spines["left"].set_color('none')
    ax1.spines["left"].set_visible(False)

    if legend:
        gslegend = mgs.GridSpecFromSubplotSpec(
            5, 1, subplot_spec=gs[2], wspace=0.0, hspace=0.0)
        epiavg = func_data.mean(3)
        epinii = nb.Nifti1Image(epiavg, img_nii.affine, img_nii.header)
        segnii = nb.Nifti1Image(lut[atlaslabels.astype(int)], epinii.affine, epinii.header)
        segnii.set_data_dtype('uint8')

        nslices = epiavg.shape[-1]
        coords = np.linspace(int(0.10 * nslices), int(0.95 * nslices), 5).astype(np.uint8)
        for i, c in enumerate(coords.tolist()):
            ax2 = plt.subplot(gslegend[i])
            plot_img(segnii, bg_img=epinii, axes=ax2, display_mode='z',
                     annotate=False, cut_coords=[c], threshold=0.1, cmap=mycolors,
                     interpolation='nearest')

    return [ax0, ax1], gs


data_dir = '/home/anibalsolon/cpac_tests/victor/outputs/o/pipeline_CUNMET_preprocessing_v2_bet/1016_ses-1'

csf_mask = data_dir + '/anatomical_csf_mask/segment_seg_0_maths.nii.gz'
csf_mask_mni = data_dir + '/anatomical_csf_mask_to_mni/segment_seg_0_maths.nii.gz'
wm_mask = data_dir + '/anatomical_wm_mask/segment_seg_2_maths.nii.gz'
wm_mask_mni = data_dir + '/anatomical_wm_mask_to_mni/segment_seg_2_maths.nii.gz'
gm_mask = data_dir + '/anatomical_gm_mask/segment_seg_1_maths.nii.gz'
gm_mask_mni = data_dir + '/anatomical_gm_mask_to_mni/segment_seg_1_maths.nii.gz'

func = data_dir + '/functional_to_standard/_scan_func-1/_compcor_ncomponents_5_selector_pc10.linear1.wm1.global0.motion1.quadratic1.gm0.compcor1.csf1/_threshold_0.3/_target_angle_deg_90/_bandpass_freqs_0.01.0.1/scrubbed_preprocessed_warp.nii.gz'
warp = data_dir + '/anatomical_to_mni_nonlinear_xfm/s1016_fieldwarp.nii.gz'

template = '/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz'

# /home/anibalsolon/cpac_tests/victor/outputs/o/pipeline_CUNMET_preprocessing_v2_bet/1016_ses-1/anatomical_to_mni_linear_xfm/s1016_calc_flirt.mat
# /home/anibalsolon/cpac_tests/victor/outputs/o/pipeline_CUNMET_preprocessing_v2_bet/1016_ses-1/functional_nuisance_regressors/_scan_func-1/_compcor_ncomponents_5_selector_pc10.linear1.wm1.global0.motion1.quadratic1.gm0.compcor1.csf1/nuisance_regressors.mat
# /home/anibalsolon/cpac_tests/victor/outputs/o/pipeline_CUNMET_preprocessing_v2_bet/1016_ses-1/functional_to_anat_linear_xfm/_scan_func-1/s1016_tshift_resample_volreg_calc_tstat_flirt.mat
# /home/anibalsolon/cpac_tests/victor/outputs/o/pipeline_CUNMET_preprocessing_v2_bet/1016_ses-1/mni_to_anatomical_linear_xfm/s1016_calc_flirt_inv.mat

import os
import subprocess

warp_cmd = 'applywarp -i {input} -r {template} -o {output} -w {warp}'

if not os.path.isfile(csf_mask_mni):
    os.makedirs(os.path.dirname(csf_mask_mni))
    subprocess.call(warp_cmd.format(
        input=csf_mask,
        template=template,
        output=csf_mask_mni,
        warp=warp
    ), shell=True)

if not os.path.isfile(wm_mask_mni):
    os.makedirs(os.path.dirname(wm_mask_mni))
    subprocess.call(warp_cmd.format(
        input=wm_mask,
        template=template,
        output=wm_mask_mni,
        warp=warp
    ), shell=True)

if not os.path.isfile(gm_mask_mni):
    os.makedirs(os.path.dirname(gm_mask_mni))
    subprocess.call(warp_cmd.format(
        input=gm_mask,
        template=template,
        output=gm_mask_mni,
        warp=warp
    ), shell=True)

func_img = nb.load(func)

gm_mask_data = (resample_to_img(nb.load(gm_mask_mni), func_img).get_data() > 0.5).astype(np.int_) * 1

wm_mask_data = (resample_to_img(nb.load(wm_mask_mni), func_img).get_data() > 0.5).astype(np.int_) * 2
wm_mask_data[gm_mask_data > 0] = 0.0

csf_mask_data = (resample_to_img(nb.load(csf_mask_mni), func_img).get_data() > 0.5).astype(np.int_) * 3
csf_mask_data[gm_mask_data > 0] = 0.0
csf_mask_data[wm_mask_data > 0] = 0.0

plot_carpet(func_img, csf_mask_data + wm_mask_data + gm_mask_data)
plt.show()