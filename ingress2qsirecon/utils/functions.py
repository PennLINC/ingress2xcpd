#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""
Convenience functions and information for ingress2qsirecon
"""

import os
import re
from pathlib import Path
from warnings import warn

import nibabel as nb

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
            "MNI_ref": ["MNINonLinear", "T1w_restore_brain.nii.gz"],
        },
        "MNI_TEMPLATE": "MNI152NLin6Asym",
        "DIR_PATTERN": re.compile(r"(\d+)"),
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
            # "MNI_ref": ["T1", "T1_unbiased_brain.nii.gz"],
        },
        "MNI_TEMPLATE": "MNI152NLin6Asym",
        "DIR_PATTERN": re.compile(r"(\d+)_(\d+)_(\d+)"),
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

    # BIDS-ify required files
    bids_dwi_file = os.path.join(bids_base, "dwi", sub_session_string + "_space-T1w_desc-preproc_dwi.nii.gz")
    bids_bval_file = os.path.join(bids_base, "dwi", sub_session_string + "_space-T1w_desc-preproc_dwi.bval")
    bids_bvec_file = os.path.join(bids_base, "dwi", sub_session_string + "_space-T1w_desc-preproc_dwi.bvec")
    bids_b_file = os.path.join(bids_base, "dwi", sub_session_string + "_space-T1w_desc-preproc_dwi.b")
    bids_bmtxt_file = os.path.join(bids_base, "dwi", sub_session_string + "_space-T1w_desc-preproc_dwi.bmtxt")
    bids_dwiref_file = os.path.join(bids_base, "dwi", sub_session_string + "_space-T1w_dwiref.nii.gz")

    bids_file_paths = {
        "bids_dwi": Path(bids_dwi_file),
        "bids_bvals": Path(bids_bval_file),
        "bids_bvecs": Path(bids_bvec_file),
        "bids_b": Path(bids_b_file),
        "bids_bmtxt": Path(bids_bmtxt_file),
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


def create_layout(input_dir: Path, output_dir: Path, input_pipeline: str, participant_label: list):
    """Find all valid directories under input_dir.
    Parameters
    ----------
    input_dir : :obj:`pathlib.Path`
        The path to the input directory.
    output_dir : :obj:`pathlib.Path`
        The path to the output directory.
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
            warn(
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
        warn(
            f"Requested participant(s) {missing_particpants} not found in layout, please confirm their data exists and are properly organized."
        )

    # Stop code if layout is empty
    if len(layout) == 0:
        raise ValueError("No subjects found in layout.")

    return layout


def to_lps(input_img, new_axcodes=("L", "P", "S")):
    if isinstance(input_img, str):
        input_img = nb.load(input_img)
    input_axcodes = nb.aff2axcodes(input_img.affine)
    # Is the input image oriented how we want?
    if not input_axcodes == new_axcodes:
        # Re-orient
        input_orientation = nb.orientations.axcodes2ornt(input_axcodes)
        desired_orientation = nb.orientations.axcodes2ornt(new_axcodes)
        transform_orientation = nb.orientations.ornt_transform(input_orientation, desired_orientation)
        reoriented_img = input_img.as_reoriented(transform_orientation)
        return reoriented_img
    else:
        return input_img
