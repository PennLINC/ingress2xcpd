#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import re
import shutil
from pathlib import Path
from textwrap import indent

import nibabel as nb
import numpy as np
from nipype import logging
from nipype.interfaces.base import (
    BaseInterface,
    BaseInterfaceInputSpec,
    Directory,
    File,
    SimpleInterface,
    Str,
    TraitedSpec,
    traits,
)

LOGGER = logging.getLogger("nipype.interface")

# Files: file paths relative to subject input folder
# MNI: MNI template version
# DIR_PATTERN: regex pattern for subject ID within input folder
PIPELINE_INFO = {
    "hcpya": {
        "files": {
            "bvals": ["T1w", "Diffusion", "bvals"],
            "bvecs": ["T1w", "Diffusion", "bvecs"],
            "dwi": ["T1w", "Diffusion", "data.nii.gz"],
            "t1w_brain": ["T1w", "T1w_acpc_dc_restore_brain.nii.gz"],
            "brain_mask": ["T1w", "brainmask_fs.nii.gz"],
            "subject2MNI": ["MNINonLinear", "xfms", "acpc_dc2standard.nii.gz"],
            "MNI2subject": ["MNINonLinear", "xfms", "standard2acpc_dc.nii.gz"],
        },
        "MNI_TEMPLATE": "MNI152NLin6Asym",
        "DIR_PATTERN": re.compile("(\d+)"),
    },
    "ukb": {
        "files": {
            "bvals": ["DTI", "dMRI", "dMRI", "bvals"],
            "bvecs": ["DTI", "dMRI", "dMRI", "bvecs"],
            "dwi": ["DTI", "dMRI", "dMRI", "data_ud.nii.gz"],
            "dwiref": ["DTI", "dMRI", "dMRI", "dti_FA.nii.gz"],
            "t1w_brain": ["T1", "T1_unbiased_brain.nii.gz"],
            "brain_mask": ["T1", "T1_brain_mask.nii.gz"],
            # TODO: Add UKB XFM path
            # "subject2MNI": ["MNINonLinear", "xfms", "acpc_dc2standard.nii.gz"], # Note this is MNI152NLin6Asym
            # "MNI2subject": ["MNINonLinear", "xfms", "standard2acpc_dc.nii.gz"], # Note this is MNI152NLin6Asym
        },
        "MNI_TEMPLATE": "MNI152NLin6Asym",
        "DIR_PATTERN": re.compile("(\d+)_(\d+)_(\d+)"),
    },
}

# List of files that are required for any recon (some pipelines require others too, this is minimum)
required_for_any_recon = [
    "bvals",
    "bvecs",
    "dwi",
]


def get_file_paths(subject_dir: Path, input_pipeline: str) -> dict:
    """Get file paths within a subject directory.

    Parameters
    ----------
    subject_dir : :obj:`pathlib.Path`
        The path to the ukb subject directory.
    input_pipeline : :obj:`str`
        The input pipeline used to create the subject directory.

    Returns
    -------
    file_paths : :obj:`dict`
        A dictionary of file paths.
    """

    # Get organization for input pipeline
    organization = PIPELINE_INFO[input_pipeline]['files']

    # Get and return file paths
    file_paths = {}
    for key, value in organization.items():
        file_paths[key] = subject_dir / Path(*value)
    return file_paths


def make_bids_file_paths(subject_layout: dict) -> dict:
    """Get file paths within a subject directory.

    Parameters
    ----------
    subject_layout : :obj:`dict`
        A dictionary of subject information from the CreateLayout function.

    Returns
    -------
    bids_file_paths : :obj:`dict`
        A dictionary of BIDS-ified file paths.
    """
    bids_base = str(subject_layout["bids_base"])
    subject = str(subject_layout["subject"])
    session = subject_layout["session"]
    if session == None:
        sub_session_string = f"sub-{subject}"
    else:
        sub_session_string = f"sub-{subject}_ses-{session}"
    mni_template = str(subject_layout["MNI_template"])

    bids_dwi_file = os.path.join(bids_base, "dwi", sub_session_string + "_dwi.nii.gz")
    bids_bval_file = os.path.join(bids_base, "dwi", sub_session_string + "_dwi.bval")
    bids_bvec_file = os.path.join(bids_base, "dwi", sub_session_string + "_dwi.bvec")
    bids_b_file = os.path.join(bids_base, "dwi", sub_session_string + "_dwi.b")
    bids_dwiref_file = os.path.join(bids_base, "dwi", sub_session_string + "_dwiref.nii.gz")

    bids_file_paths = {
        "bids_dwi": Path(bids_dwi_file),
        "bids_bvals": Path(bids_bval_file),
        "bids_bvecs": Path(bids_bvec_file),
        "bids_b": Path(bids_b_file),
        "bids_dwiref": Path(bids_dwiref_file),
    }

    # Now for optional files
    if subject_layout['t1w_brain']:
        bids_t1w_brain = os.path.join(bids_base, "anat", sub_session_string + "_desc-preproc_T1w.nii.gz")
        bids_file_paths.update({"bids_t1w_brain": Path(bids_t1w_brain)})
    if subject_layout['brain_mask']:
        bids_brain_mask = os.path.join(bids_base, "anat", sub_session_string + "_desc-brain_mask.nii.gz")
        bids_file_paths.update({"bids_brain_mask": Path(bids_brain_mask)})
    if subject_layout['subject2MNI']:
        bids_subject2MNI = os.path.join(
            bids_base, "anat", sub_session_string + f"_from-T1w_to-{mni_template}_mode-image_xfm.h5"
        )
        bids_file_paths.update({"bids_subject2MNI": Path(bids_subject2MNI)})
    if subject_layout['MNI2subject']:
        bids_MNI2subject = os.path.join(
            bids_base, "anat", sub_session_string + f"_from-{mni_template}_to-T1w_mode-image_xfm.h5"
        )
        bids_file_paths.update({"bids_MNI2subject": Path(bids_MNI2subject)})

    return bids_file_paths


class _CreateLayoutInputSpec(BaseInterfaceInputSpec):
    input_dir = Directory(exists=True, mandatory=True, desc="Path to the input directory")
    output_dir = Directory(exists=True, mandatory=True, desc="Path to the output directory")
    input_pipeline = Str(mandatory=True, desc="Name of the input pipeline")
    participant_label = traits.List(Str, desc="List of participant labels to search for")


class _CreateLayoutOutputSpec(TraitedSpec):
    layout = traits.List(traits.Dict, desc="List of dictionaries with subject ID, session ID, and paths")


class CreateLayout(SimpleInterface):
    input_spec = _CreateLayoutInputSpec
    output_spec = _CreateLayoutOutputSpec

    def _run_interface(self, runtime):
        input_dir = Path(self.inputs.input_dir)
        output_dir = Path(self.inputs.output_dir)
        input_pipeline = self.inputs.input_pipeline
        participant_label = self.inputs.participant_label
        pattern = PIPELINE_INFO[input_pipeline]["DIR_PATTERN"]  # search pattern for subject ID, session etc
        MNI_template = PIPELINE_INFO[input_pipeline]["MNI_TEMPLATE"]

        layout = []
        for potential_dir in input_dir.iterdir():  # loop through all directories in input_dir
            if participant_label and not potential_dir.name.startswith(tuple(participant_label)):
                # Skip if folder in loop does not start with an expected participant label
                continue

            match = re.match(pattern, potential_dir.name)
            if not match:
                # Skip if subject folder does not match expected pattern
                continue

            # If passes all checks, add to layout
            if input_pipeline == "hcpya":
                subject = match.group(1)
                renamed_ses = None
            elif input_pipeline == "ukb":
                subject, ses_major, ses_minor = match.groups()
                renamed_ses = "%02d%02d" % (int(ses_major), int(ses_minor))

            # Make BIDS base organization
            bids_base = output_dir / f"sub-{subject}"
            if renamed_ses:
                bids_base = bids_base / f"ses-{renamed_ses}"

            file_paths = get_file_paths(potential_dir, input_pipeline)
            # check if any required files do not exist
            missing_for_any_recon = [
                file_type for file_type in required_for_any_recon if not os.path.isfile(file_paths[file_type])
            ]
            if missing_for_any_recon:
                LOGGER.warning(
                    f"Required files missing for any recon: {missing_for_any_recon}. "
                    "These are expected at the following locations: "
                    f"{[file_type + ': ' + str(file_paths[file_type]) for file_type in missing_for_any_recon]}. "
                    f"Skipping subject {subject}."
                )
                continue

            subject_layout = {
                "subject": subject,
                "session": renamed_ses,
                "path": Path(potential_dir),
                "bids_base": bids_base,
                "MNI_template": MNI_template,
            }
            # Add file paths to subject layout if they exist
            subject_layout.update({file_type: path for file_type, path in file_paths.items() if os.path.exists(path)})
            # Make BIDS-file path names
            subject_layout.update(make_bids_file_paths(subject_layout))
            # Save out layout
            layout.append(subject_layout)

        # Sort layout by subject ID
        layout = sorted(layout, key=lambda x: x["subject"])

        # Raise warnings for requested subjects not in layout
        missing_particpants = sorted(set(participant_label) - set([subject["subject"] for subject in layout]))
        if missing_particpants:
            LOGGER.warning(
                f"Requested participant(s) {missing_particpants} not found in layout, please confirm their data exists and are properly organized."
            )

        # Stop code if layout is empty
        if len(layout) == 0:
            raise ValueError("No subjects found in layout.")

        self._results = {}
        self._results['layout'] = layout
        return runtime


class _ValidateImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="input image")
    out_file = File(exists=False, desc="validated image")


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
        out_report = os.path.join(runtime.cwd, "report.html")

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

        validator = ValidateImage(in_file=dwi_in_file, out_file=dwi_out_file)
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
            self._results["dwi_file"] = dwi_out_file

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
                self._results["bvec_file"] = bvec_out_file

        else:
            LOGGER.info("Not applying reorientation to %s: already in %s", dwi_in_file, orientation)
            self._results["dwi_file"] = dwi_out_file
            # Copy and rename bvecs
            shutil.copy(bvec_in_file, bvec_out_file)
            self._results["bvec_file"] = bvec_out_file

        # Copy and rename bvals
        shutil.copy(bval_in_file, bval_out_file)
        self._results["bval_file"] = bval_out_file

        return runtime


####################################


class _DirnameToBidsInputSpec(BaseInterfaceInputSpec):
    subject_dir = File(exists=True, mandatory=True, desc="Path to the subject directory")
    input_pipeline = Str(mandatory=True, desc="Name of the input pipeline (e.g. 'hcpya', 'ukb')")


class _DirnameToBidsOutputSpec(TraitedSpec):
    bids_id = Str(desc="BIDS subject ID and session ID (if applicable)")


class DirnameToBids(BaseInterface):
    input_spec = _DirnameToBidsInputSpec
    output_spec = _DirnameToBidsOutputSpec

    def _run_interface(self, runtime):
        subject_dir = Path(self.inputs.subject_dir)
        input_pipeline = self.inputs.input_pipeline

        pattern = DIR_PATTERNS[input_pipeline]
        match = re.match(pattern, subject_dir.name)

        if input_pipeline == "hcpya":
            subject = match.group(1)
            bids_id = f"sub-{subject}"
        elif input_pipeline == "ukb":
            subject, ses_major, ses_minor = match.groups()
            renamed_ses = "%02d%02d" % (int(ses_major), int(ses_minor))
            bids_id = f"sub-{subject}_ses-{renamed_ses}"

        self._results['bids_id'] = bids_id
        return runtime


class _GetFilePathsInputSpec(BaseInterfaceInputSpec):
    subject_dir = Directory(exists=True, mandatory=True, desc="Path to the subject directory")
    input_pipeline = Str(mandatory=True, desc="Input pipeline used to create the subject directory")
    subject_label = Str(mandatory=True, desc="Subject folder identifier")


class _GetFilePathsOutputSpec(TraitedSpec):
    file_paths = traits.Dict(desc="Dictionary of file paths")


class GetFilePaths(BaseInterface):
    input_spec = _GetFilePathsInputSpec
    output_spec = _GetFilePathsOutputSpec

    def _run_interface(self, runtime):
        ORGANIZATIONS = {
            "hcpya": {
                "bvals": ["T1w", "Diffusion", "bvals"],
                "bvecs": ["T1w", "Diffusion", "bvecs"],
                "dwi": ["T1w", "Diffusion", "data.nii.gz"],
                "t1w_brain": ["T1w", "T1w_acpc_dc_restore_brain.nii.gz"],
                "brain_mask": ["T1w", "brainmask_fs.nii.gz"],
                "subject2MNI": ["MNINonLinear", "xfms", "acpc_dc2standard.nii.gz"],
                "freesurfer": ["T1w", self.inputs.subject_label],
            },
            "ukb": {
                "bvals": ["DTI", "dMRI", "dMRI", "bvals"],
                "bvecs": ["DTI", "dMRI", "dMRI", "bvecs"],
                "dwi": ["DTI", "dMRI", "dMRI", "data_ud.nii.gz"],
                "dwiref": ["DTI", "dMRI", "dMRI", "dti_FA.nii.gz"],
                # Additional paths can be added here
            },
        }

        organization = ORGANIZATIONS[self.inputs.input_pipeline]
        subject_dir = Path(self.inputs.subject_dir)
        file_paths = {key: subject_dir / Path(*value) for key, value in organization.items()}

        self._results['file_paths'] = file_paths
        return runtime


class _CollectParticipantsInputSpec(BaseInterfaceInputSpec):
    layout = traits.List(traits.Dict, mandatory=True, desc="List of dictionaries with subject info")
    participant_label = traits.List(Str, desc="List of participant labels to search for")


class _CollectParticipantsOutputSpec(TraitedSpec):
    participants = traits.List(Str, desc="List of valid participant labels")


class CollectParticipants(SimpleInterface):
    input_spec = _CollectParticipantsInputSpec
    output_spec = _CollectParticipantsOutputSpec

    def _run_interface(self, runtime):
        layout = self.inputs.layout
        participant_label = self.inputs.participant_label

        all_participants = set([spec["subject"] for spec in layout])

        # if not participant_label:
        if len(participant_label) == 0:
            participants = sorted(all_participants)
        else:
            participant_label = set(participant_label)

            found_labels = sorted(participant_label & all_participants)
            requested_but_missing = sorted(participant_label - all_participants)

            if requested_but_missing:
                raise Exception(
                    f"Requested subjects [{', '.join(requested_but_missing)}] do not have complete directories."
                )
            if not found_labels:
                raise Exception("No complete directories were found")

            participants = sorted(found_labels)

        self._results = {}
        self._results['participants'] = participants
        return runtime


#####################################


def dirname_to_bids(subject_dir: Path, input_pipeline: str) -> str:
    """Convert a subject directory name to a BIDS subject ID.

    Parameters
    ----------
    subject_dir : :obj:`pathlib.Path`
        The path to the subject directory.
    input_pipeline : :obj:`str`
        The name of the input pipeline (e.g. 'hcpya', 'ukb')

    Returns
    -------
    bids_id : :obj:`str`
        The BIDS subject ID and session ID (if applicable), as a string.
    """

    pattern = DIR_PATTERNS[input_pipeline]
    match = re.match(pattern, subject_dir.name)
    if input_pipeline == "hcpya":
        subject = match.group(1)
        return f"sub-{subject}"
    elif input_pipeline == "ukb":
        subject, ses_major, ses_minor = match.groups()
        renamed_ses = "%02d%02d" % (int(ses_major), int(ses_minor))
        return f"sub-{subject}_ses-{renamed_ses}"


def create_layout(input_dir, input_pipeline, participant_label=None):
    """Find all valid directories under input_dir.
    Parameters
    ----------
    input_dir : :obj:`pathlib.Path`
        The path to the input directory.
    input_pipeline : :obj:`str`
        The name of the input pipeline (e.g. 'hcpya', 'ukb')
    participant_label : :obj:`list` of :obj:`str`
        A list of participant labels to search for.
    Returns
    -------
    layout : :obj:`list` of :obj:`dict`
        A list of dictionaries containing the subject ID, session ID (if applicable), path to the subject directory,
        and the path to the fake dwi file.
    """
    # find directories starting with a number. These are the candidate directories
    pattern = DIR_PATTERNS[input_pipeline]
    layout = []
    for potential_dir in Path(input_dir).iterdir():
        if participant_label and not potential_dir.name.startswith(participant_label):
            continue

        match = re.match(pattern, potential_dir.name)
        if not match:
            continue

        if missing_from_subject_dir(potential_dir):
            continue

        subject = match.groups()[0]
        fake_dwi_file = f"/bids/sub-{subject}/dwi/sub-{subject}_dwi.nii.gz"
        layout.append(
            {
                "subject": subject,
                "path": potential_dir,
                "bids_dwi_file": fake_dwi_file,
            }
        )

    return layout


def missing_from_subject_dir(subject_dir: Path, input_pipeline: str, subject_label: str) -> list:
    """Check for missing files in a subject directory.

    Parameters
    ----------
    subject_dir : :obj:`pathlib.Path`
        The path to the ukb subject directory.
    input_pipeline : :obj:`str`
        The name of the input pipeline (e.g. 'hcpya', 'ukb')
    subject_label : :obj:`str`
        The subject folder identifier.

    Returns
    -------
    missing_files : :obj:`list` of :obj:`str`
        A list of missing files.
    """
    # Files needed to run recon workflows
    # Anatomical information not necessary for some workflows
    fields = [
        "bvals",
        "bvecs",
        "dwi",
        "dwiref",
        "t1w_brain",
        "brain_mask",
        "subject2MNI",
        "freesurfer",
    ]
    required_for_any_recon = [
        "bvals",
        "bvecs",
        "dwi",
    ]
    file_paths = get_file_paths(subject_dir, input_pipeline, subject_label)

    pass
    return None


def collect_participants(layout, participant_label):
    """Collect all valid participants.

    Parameters
    ----------
    layout : :obj:`list` of :obj:`dict`
        A list of dictionaries containing the subject ID, session ID (if applicable), path to the directory,
        and the path to the fake dwi file.
    participant_label : :obj:`list` of :obj:`str`
        A list of participant labels to search for.

    Returns
    -------
    participants : :obj:`list` of :obj:`str`
        A list of participant labels.
    """
    all_participants = set([spec["subject"] for spec in layout])

    # No --participant-label was set, return all
    if not participant_label:
        return sorted(all_participants)

    if isinstance(participant_label, str):
        participant_label = [participant_label]

    participant_label = set(participant_label)

    # Remove labels not found
    found_labels = sorted(participant_label & all_participants)
    requested_but_missing = sorted(participant_label - all_participants)

    if requested_but_missing:
        raise Exception(
            "Requested subjects [{}] do not have complete directories under".format(", ".join(requested_but_missing))
        )
    if not found_labels:
        raise Exception("No complete directories were found")

    return sorted(found_labels)
