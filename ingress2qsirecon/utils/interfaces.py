"""
Nipype Interfaces for Ingress2Qsirecon
"""

import os
import shutil
from textwrap import indent

import nibabel as nb
import numpy as np
import pandas as pd
import SimpleITK as sitk
from nilearn import image as nim
from nipype import logging
from nipype.interfaces import ants
from nipype.interfaces.base import (
    BaseInterface,
    BaseInterfaceInputSpec,
    CommandLineInputSpec,
    File,
    InputMultiPath,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.interfaces.mixins import reporting
from nipype.interfaces.workbench.base import WBCommand
from niworkflows.interfaces.norm import (
    SpatialNormalization,
    _SpatialNormalizationInputSpec,
)
from niworkflows.interfaces.reportlets.base import (
    RegistrationRC,
    _SVGReportCapableInputSpec,
)

from ingress2qsirecon.utils.functions import to_lps

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


class _ConformInputSpec(BaseInterfaceInputSpec):
    in_file = File(mandatory=True, desc="Input image")
    out_file = File(mandatory=True, genfile=True, desc="Conformed image")
    target_zooms = traits.Tuple(traits.Float, traits.Float, traits.Float, desc="Target zoom information")
    target_shape = traits.Tuple(traits.Int, traits.Int, traits.Int, desc="Target shape information")
    deoblique_header = traits.Bool(False, usedfault=True)


class _ConformOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Conformed image")
    # transform = File(exists=True, desc="Conformation transform")
    # report = File(exists=True, desc='reportlet about orientation')


class Conform(SimpleInterface):
    """Conform a series of T1w images to enable merging.

    Performs two basic functions:

    1. Orient to LPS (right-left, anterior-posterior, inferior-superior)
    2. Resample to target zooms (voxel sizes) and shape (number of voxels)

    """

    input_spec = _ConformInputSpec
    output_spec = _ConformOutputSpec

    def _run_interface(self, runtime):
        # Load image, orient as LPS
        fname = self.inputs.in_file
        orig_img = nb.load(fname)
        reoriented = to_lps(orig_img)

        # Set target shape information
        target_zooms = np.array(self.inputs.target_zooms)
        target_shape = np.array(self.inputs.target_shape)
        target_span = target_shape * target_zooms

        zooms = np.array(reoriented.header.get_zooms()[:3])
        shape = np.array(reoriented.shape[:3])

        # Reconstruct transform from orig to reoriented image
        ornt_xfm = nb.orientations.inv_ornt_aff(nb.io_orientation(reoriented.affine), orig_img.shape)
        # Identity unless proven otherwise
        target_affine = reoriented.affine.copy()
        conform_xfm = np.eye(4)
        # conform_xfm = np.diag([-1, -1, 1, 1])

        xyz_unit = reoriented.header.get_xyzt_units()[0]
        if xyz_unit == "unknown":
            # Common assumption; if we're wrong, unlikely to be the only thing that breaks
            xyz_unit = "mm"

        # Set a 0.05mm threshold to performing rescaling
        atol = {"meter": 1e-5, "mm": 0.01, "micron": 10}[xyz_unit]

        # Rescale => change zooms
        # Resize => update image dimensions
        rescale = not np.allclose(zooms, target_zooms, atol=atol)
        resize = not np.all(shape == target_shape)
        if rescale or resize:
            if rescale:
                scale_factor = target_zooms / zooms
                target_affine[:3, :3] = reoriented.affine[:3, :3].dot(np.diag(scale_factor))

            if resize:
                # The shift is applied after scaling.
                # Use a proportional shift to maintain relative position in dataset
                size_factor = target_span / (zooms * shape)
                # Use integer shifts to avoid unnecessary interpolation
                offset = reoriented.affine[:3, 3] * size_factor - reoriented.affine[:3, 3]
                target_affine[:3, 3] = reoriented.affine[:3, 3] + offset.astype(int)

            data = nim.resample_img(reoriented, target_affine, target_shape).get_fdata()
            conform_xfm = np.linalg.inv(reoriented.affine).dot(target_affine)
            reoriented = reoriented.__class__(data, target_affine, reoriented.header)

        if self.inputs.deoblique_header:
            is_oblique = np.any(np.abs(nb.affines.obliquity(reoriented.affine)) > 0)
            if is_oblique:
                LOGGER.warning("Removing obliquity from image affine")
                new_affine = reoriented.affine.copy()
                new_affine[:, :-1] = 0
                new_affine[(0, 1, 2), (0, 1, 2)] = reoriented.header.get_zooms()[:3] * np.sign(
                    reoriented.affine[(0, 1, 2), (0, 1, 2)]
                )
                reoriented = nb.Nifti1Image(reoriented.get_fdata(), new_affine, reoriented.header)

        # Save image
        out_name = self.inputs.out_file
        reoriented.to_filename(out_name)

        # Image may be reoriented, rescaled, and/or resized
        if reoriented is not orig_img:

            transform = ornt_xfm.dot(conform_xfm)
            if not np.allclose(orig_img.affine.dot(transform), target_affine):
                LOGGER.warning("Check alignment of anatomical image.")

        else:
            transform = np.eye(4)

        # mat_name = fname_presuffix(fname, suffix=".mat", newpath=runtime.cwd, use_ext=False)
        # np.savetxt(mat_name, transform, fmt="%.08f")
        # self._results["transform"] = mat_name
        self._results["out_file"] = out_name

        return runtime


# Define the input specification
class _ComposeTransformsInputSpec(BaseInterfaceInputSpec):
    warp_files = InputMultiPath(File(exists=True), desc="List of warp files in .h5 format", mandatory=True)
    output_warp = File(mandatory=True, genfile=True, desc="Output composed warp file")


# Define the output specification
class _ComposeTransformsOutputSpec(TraitedSpec):
    output_warp = File(exists=True, desc="Output composed warp file")


# Define the custom Nipype interface
class ComposeTransforms(BaseInterface):

    input_spec = _ComposeTransformsInputSpec
    output_spec = _ComposeTransformsOutputSpec

    def _run_interface(self, runtime):
        # Create a CompositeTransform object
        composite_transform = sitk.CompositeTransform(3)

        # Iterate over the list of warp files and add them to the composite transform
        for warp_file in self.inputs.warp_files:
            transform = sitk.ReadTransform(warp_file)
            try:
                # If composite, add each transform in the list
                for i in range(transform.GetNumberOfTransforms()):
                    composite_transform.AddTransform(transform.GetNthTransform(i))
            except:
                # If not, just add the transform
                composite_transform.AddTransform(transform)

        # Write the composite transform to the temporary file
        sitk.WriteTransform(composite_transform, self.inputs.output_warp)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["output_warp"] = os.path.abspath(self.inputs.output_warp)
        return outputs


class _FSLBVecsToTORTOISEBmatrixInputSpec(BaseInterfaceInputSpec):
    bvals_file = File(exists=True, desc='Full path to bvals file', mandatory=True)
    bvecs_file = File(exists=True, desc='Full path to bvecs file', mandatory=True)
    bmtxt_file = File(mandatory=True, desc='Output B-matrix file', genfile=True)


class _FSLBVecsToTORTOISEBmatrixOutputSpec(TraitedSpec):
    bmtxt_file = File(exists=True, genfile=True, desc='Output B-matrix file')


class FSLBVecsToTORTOISEBmatrix(BaseInterface):
    input_spec = _FSLBVecsToTORTOISEBmatrixInputSpec
    output_spec = _FSLBVecsToTORTOISEBmatrixOutputSpec

    def _run_interface(self, runtime):
        bvals_file = self.inputs.bvals_file
        bvecs_file = self.inputs.bvecs_file

        # Load bvals and bvecs
        try:
            bvals = np.loadtxt(bvals_file)
        except OSError:
            raise RuntimeError(f"Bvals file does not exist: {bvals_file}")

        try:
            bvecs = np.loadtxt(bvecs_file)
        except OSError:
            raise RuntimeError(f"Bvecs file does not exist: {bvecs_file}")

        # Ensure bvecs has 3 rows and bvals has 1 row
        if bvecs.shape[0] != 3:
            bvecs = bvecs.T
        if bvals.shape[0] != 1:
            bvals = bvals.reshape(1, -1)

        Nvols = bvecs.shape[1]
        Bmatrix = np.zeros((Nvols, 6))

        for i in range(Nvols):
            vec = bvecs[:, i].reshape(3, 1)
            nrm = np.linalg.norm(vec)
            if nrm > 1e-3:
                vec /= nrm

            mat = bvals[0, i] * np.dot(vec, vec.T)
            Bmatrix[i, 0] = mat[0, 0]
            Bmatrix[i, 1] = 2 * mat[0, 1]
            Bmatrix[i, 2] = 2 * mat[0, 2]
            Bmatrix[i, 3] = mat[1, 1]
            Bmatrix[i, 4] = 2 * mat[1, 2]
            Bmatrix[i, 5] = mat[2, 2]

        dirname = os.path.dirname(bvals_file)
        if dirname == "":
            dirname = "."

        np.savetxt(self.inputs.bmtxt_file, Bmatrix)

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['bmtxt_file'] = self.inputs.bmtxt_file
        return outputs


class _MRTrixGradientTableInputSpec(BaseInterfaceInputSpec):
    bval_file = File(exists=True, mandatory=True)
    bvec_file = File(exists=True, mandatory=True)
    b_file_out = File(genfile=True, mandatory=True)


class _MRTrixGradientTableOutputSpec(TraitedSpec):
    b_file_out = File(exists=True)


class MRTrixGradientTable(BaseInterface):
    input_spec = _MRTrixGradientTableInputSpec
    output_spec = _MRTrixGradientTableOutputSpec

    def _run_interface(self, runtime):
        _convert_fsl_to_mrtrix(self.inputs.bval_file, self.inputs.bvec_file, self.inputs.b_file_out)
        return runtime

    def _list_outputs(self):
        outputs = self.output_spec().get()
        # Return the output filename which Nipype generates
        outputs['b_file_out'] = self.inputs.b_file_out
        return outputs


def _convert_fsl_to_mrtrix(bval_file, bvec_file, output_fname):
    vecs = np.loadtxt(bvec_file)
    vals = np.loadtxt(bval_file)
    gtab = np.column_stack([vecs.T, vals]) * np.array([-1, -1, 1, 1])
    np.savetxt(output_fname, gtab, fmt=["%.8f", "%.8f", "%.8f", "%d"])


class RobustMNINormalizationInputSpecRPT(
    _SVGReportCapableInputSpec,
    _SpatialNormalizationInputSpec,
):
    # Template orientation.
    orientation = traits.Enum(
        "LPS",
        mandatory=True,
        usedefault=True,
        desc="modify template orientation (should match input image)",
    )


class RobustMNINormalizationOutputSpecRPT(
    reporting.ReportCapableOutputSpec,
    ants.registration.RegistrationOutputSpec,
):
    # Try to work around TraitError of "undefined 'reference_image' attribute"
    reference_image = traits.File(desc="the output reference image")


class RobustMNINormalizationRPT(RegistrationRC, SpatialNormalization):
    input_spec = RobustMNINormalizationInputSpecRPT
    output_spec = RobustMNINormalizationOutputSpecRPT

    def _post_run_hook(self, runtime):
        # We need to dig into the internal ants.Registration interface
        self._fixed_image = self._get_ants_args()["fixed_image"]
        if isinstance(self._fixed_image, (list, tuple)):
            self._fixed_image = self._fixed_image[0]  # get first item if list

        if self._get_ants_args().get("fixed_image_mask") is not None:
            self._fixed_image_mask = self._get_ants_args().get("fixed_image_mask")
        self._moving_image = self.aggregate_outputs(runtime=runtime).warped_image
        LOGGER.info(
            "Report - setting fixed (%s) and moving (%s) images",
            self._fixed_image,
            self._moving_image,
        )

        return super(RobustMNINormalizationRPT, self)._post_run_hook(runtime)
