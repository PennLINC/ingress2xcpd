"""
Nipype Interfaces for Ingress2Qsirecon
"""

import os
import shutil
from textwrap import indent

import nibabel as nb
import numpy as np
import SimpleITK as sitk
from nilearn import image as nim
from nipype import logging
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    CommandLineInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.interfaces.workbench.base import WBCommand

LOGGER = logging.getLogger("nipype.interface")


class _ValidateImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="input image")
    out_file = File(mandatory=True, desc="validated image", genfile=True)
    out_report = File(mandatory=True, desc="HTML segment containing warning", genfile=True)


class _ValidateImageOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="validated image")
    out_report = File(exists=True, desc="HTML segment containing warning")


class ValidateImage(SimpleInterface):
    """
    Check the correctness of x-form headers (matrix and code)
    This interface implements the `following logic
    <https://github.com/poldracklab/fmriprep/issues/873#issuecomment-349394544>`_:
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    | valid quaternions | `qform_code > 0` | `sform_code > 0` | `qform == sform` \
| actions                                        |
    +===================+==================+==================+==================\
+================================================+
    | True              | True             | True             | True             \
| None                                           |
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    | True              | True             | False            | *                \
| sform, scode <- qform, qcode                   |
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    | *                 | *                | True             | False            \
| qform, qcode <- sform, scode                   |
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    | *                 | False            | True             | *                \
| qform, qcode <- sform, scode                   |
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    | *                 | False            | False            | *                \
| sform, qform <- best affine; scode, qcode <- 1 |
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    | False             | *                | False            | *                \
| sform, qform <- best affine; scode, qcode <- 1 |
    +-------------------+------------------+------------------+------------------\
+------------------------------------------------+
    """

    input_spec = _ValidateImageInputSpec
    output_spec = _ValidateImageOutputSpec

    def _run_interface(self, runtime):
        img = nb.load(self.inputs.in_file)
        out_report = self.inputs.out_report

        # Retrieve xform codes
        sform_code = int(img.header._structarr["sform_code"])
        qform_code = int(img.header._structarr["qform_code"])

        # Check qform is valid
        valid_qform = False
        try:
            qform = img.get_qform()
            valid_qform = True
        except ValueError:
            pass

        sform = img.get_sform()
        if np.linalg.det(sform) == 0:
            valid_sform = False
        else:
            RZS = sform[:3, :3]
            zooms = np.sqrt(np.sum(RZS * RZS, axis=0))
            valid_sform = np.allclose(zooms, img.header.get_zooms()[:3])

        # Matching affines
        matching_affines = valid_qform and np.allclose(qform, sform)

        # Both match, qform valid (implicit with match), codes okay -> do nothing, empty report
        if matching_affines and qform_code > 0 and sform_code > 0:
            self._results["out_file"] = self.inputs.in_file
            open(out_report, "w").close()
            self._results["out_report"] = out_report
            return runtime

        # A new file will be written
        out_fname = self.inputs.out_file
        self._results["out_file"] = out_fname

        # Row 2:
        if valid_qform and qform_code > 0 and (sform_code == 0 or not valid_sform):
            img.set_sform(qform, qform_code)
            warning_txt = "Note on orientation: sform matrix set"
            description = """\
<p class="elem-desc">The sform has been copied from qform.</p>
"""
        # Rows 3-4:
        # Note: if qform is not valid, matching_affines is False
        elif (valid_sform and sform_code > 0) and (not matching_affines or qform_code == 0):
            img.set_qform(img.get_sform(), sform_code)
            warning_txt = "Note on orientation: qform matrix overwritten"
            description = """\
<p class="elem-desc">The qform has been copied from sform.</p>
"""
            if not valid_qform and qform_code > 0:
                warning_txt = "WARNING - Invalid qform information"
                description = """\
<p class="elem-desc">
    The qform matrix found in the file header is invalid.
    The qform has been copied from sform.
    Checking the original qform information from the data produced
    by the scanner is advised.
</p>
"""
        # Rows 5-6:
        else:
            affine = img.header.get_base_affine()
            img.set_sform(affine, nb.nifti1.xform_codes["scanner"])
            img.set_qform(affine, nb.nifti1.xform_codes["scanner"])
            warning_txt = "WARNING - Missing orientation information"
            description = """\
<p class="elem-desc">
    QSIRecon could not retrieve orientation information from the image header.
    The qform and sform matrices have been set to a default, LAS-oriented affine.
    Analyses of this dataset MAY BE INVALID.
</p>
"""
        snippet = '<h3 class="elem-title">%s</h3>\n%s:\n\t %s\n' % (
            warning_txt,
            self.inputs.in_file,
            description,
        )
        # Store new file and report
        img.to_filename(out_fname)
        with open(out_report, "w") as fobj:
            fobj.write(indent(snippet, "\t" * 3))

        self._results["out_report"] = out_report
        return runtime


class _ConformDwiInputSpec(BaseInterfaceInputSpec):
    dwi_in_file = File(mandatory=True, desc="dwi image")
    bval_in_file = File(exists=True)
    bvec_in_file = File(exists=True)
    dwi_out_file = File(desc="conformed dwi image")
    bval_out_file = File(desc="conformed bval file")
    bvec_out_file = File(desc="conformed bvec file")
    orientation = traits.Enum("LPS", "LAS", default="LPS", usedefault=True)


class _ConformDwiOutputSpec(TraitedSpec):
    dwi_out_file = File(exists=True, desc="conformed dwi image")
    bvec_out_file = File(exists=True, desc="conformed bvec file")
    bval_out_file = File(exists=True, desc="conformed bval file")
    out_report = File(exists=True, desc="HTML segment containing warning")


class ConformDwi(SimpleInterface):
    """Conform a series of dwi images to enable merging.
    Performs three basic functions:
    #. Orient image to requested orientation
    #. Validate the qform and sform, set qform code to 1
    #. Flip bvecs accordingly
    #. Do nothing to the bvals
    Note: This is not as nuanced as fmriprep's version
    """

    input_spec = _ConformDwiInputSpec
    output_spec = _ConformDwiOutputSpec

    def _run_interface(self, runtime):
        dwi_in_file = self.inputs.dwi_in_file
        bval_in_file = self.inputs.bval_in_file
        bvec_in_file = self.inputs.bvec_in_file
        dwi_out_file = self.inputs.dwi_out_file
        bval_out_file = self.inputs.bval_out_file
        bvec_out_file = self.inputs.bvec_out_file
        orientation = self.inputs.orientation

        validator = ValidateImage(in_file=dwi_in_file, out_file=dwi_out_file, out_report=os.getcwd() + "/report.txt")
        validated = validator.run()
        self._results["out_report"] = validated.outputs.out_report
        input_img = nb.load(validated.outputs.out_file)

        input_axcodes = nb.aff2axcodes(input_img.affine)
        # Is the input image oriented how we want?
        new_axcodes = tuple(orientation)

        if not input_axcodes == new_axcodes:
            # Re-orient
            LOGGER.info("Re-orienting %s to %s", dwi_in_file, orientation)
            input_orientation = nb.orientations.axcodes2ornt(input_axcodes)
            desired_orientation = nb.orientations.axcodes2ornt(new_axcodes)
            transform_orientation = nb.orientations.ornt_transform(input_orientation, desired_orientation)
            reoriented_img = input_img.as_reoriented(transform_orientation)
            reoriented_img.to_filename(dwi_out_file)
            self._results["dwi_out_file"] = dwi_out_file

            # Flip the bvecs
            if os.path.exists(bvec_in_file):
                LOGGER.info("Reorienting %s to %s", bvec_in_file, orientation)
                bvec_array = np.loadtxt(bvec_in_file)
                if not bvec_array.shape[0] == transform_orientation.shape[0]:
                    raise ValueError("Unrecognized bvec format")
                output_array = np.zeros_like(bvec_array)
                for this_axnum, (axnum, flip) in enumerate(transform_orientation):
                    output_array[this_axnum] = bvec_array[int(axnum)] * flip
                np.savetxt(bvec_out_file, output_array, fmt="%.8f ")
                self._results["bvec_out_file"] = bvec_out_file

        else:
            LOGGER.info("Not applying reorientation to %s: already in %s", dwi_in_file, orientation)
            self._results["dwi_out_file"] = dwi_out_file
            # Copy and rename bvecs
            if not os.path.exists(bvec_out_file):
                shutil.copy(bvec_in_file, bvec_out_file)
            self._results["bvec_out_file"] = bvec_out_file

        # Copy and rename bvals
        if not os.path.exists(bval_out_file):
            shutil.copy(bval_in_file, bval_out_file)
        self._results["bval_out_file"] = bval_out_file

        return runtime


class _ConvertWarpfieldInputSpec(CommandLineInputSpec):
    fnirt_in_xfm = File(
        exists=True,
        mandatory=True,
        argstr="-from-fnirt %s",
        position=0,
        desc="The input FNIRT warp",
    )
    fnirt_ref_file = File(
        exists=True,
        mandatory=True,
        argstr="%s",
        position=1,
        desc="The reference imag used for FNIRT",
    )
    itk_out_xfm = File(
        genfile=True,
        mandatory=True,
        argstr="-to-itk %s",
        position=2,
        desc="The output ITK warp",
    )


class _ConvertWarpfieldOutputSpec(TraitedSpec):
    itk_out_xfm = File(exists=True, desc="output CIFTI file")


class ConvertWarpfield(WBCommand):
    """
    Use the wb_command to convert a FNIRT oriented .nii.gz to an ITK .nii.gz
    """

    input_spec = _ConvertWarpfieldInputSpec
    output_spec = _ConvertWarpfieldOutputSpec
    _cmd = "wb_command -convert-warpfield"

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs["itk_out_xfm"] = os.path.abspath(self.inputs.itk_out_xfm)
        return outputs


class _NIFTItoH5InputSpec(TraitedSpec):
    xfm_nifti_in = File(exists=True, mandatory=True, desc="ITK NIFTI xfm iput")
    xfm_h5_out = File(mandatory=True, desc="ITK H5 xfm output", genfile=True)


class _NIFTItoH5OutputSpec(TraitedSpec):
    xfm_h5_out = File(exists=True, desc="output image")


class NIFTItoH5(SimpleInterface):

    input_spec = _NIFTItoH5InputSpec
    output_spec = _NIFTItoH5OutputSpec

    def _run_interface(self, runtime):
        displacement_image = sitk.ReadImage(self.inputs.xfm_nifti_in, sitk.sitkVectorFloat64, imageIO="NiftiImageIO")
        tx = sitk.DisplacementFieldTransform(displacement_image)
        sitk.WriteTransform(tx, self.inputs.xfm_h5_out)
        self._results["xfm_h5_out"] = self.inputs.xfm_h5_out

        return runtime


class _ExtractB0sInputSpec(BaseInterfaceInputSpec):
    b0_indices = traits.List()
    bval_file = File(exists=True)
    b0_threshold = traits.Int(50, usedefault=True)
    dwi_series = File(exists=True, mandatory=True)
    b0_average = File(mandatory=True, genfile=True)


class _ExtractB0sOutputSpec(TraitedSpec):
    b0_average = File(exists=True)


class ExtractB0s(SimpleInterface):
    """Extract a b0 series and a mean b0 from a dwi series."""

    input_spec = _ExtractB0sInputSpec
    output_spec = _ExtractB0sOutputSpec

    def _run_interface(self, runtime):
        output_mean_fname = self.inputs.b0_average
        bvals = np.loadtxt(self.inputs.bval_file)
        indices = np.flatnonzero(bvals < self.inputs.b0_threshold)
        if indices.size == 0:
            raise ValueError("No b<%d images found" % self.inputs.b0_threshold)

        new_data = nim.index_img(self.inputs.dwi_series, indices)
        if new_data.ndim > 3:
            mean_image = nim.math_img("img.mean(3)", img=new_data)
            mean_image.to_filename(output_mean_fname)
        else:
            new_data.to_filename(output_mean_fname)

        self._results["b0_average"] = output_mean_fname

        return runtime
