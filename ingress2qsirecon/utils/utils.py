#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import re
from pathlib import Path

DIR_PATTERNS = {
        'hcpya': re.compile("(\d+)"), # HCP is just a number
        'ukb': re.compile("(\d+)_(\d+)_(\d+)"), # UKB has three numbers for subject, session major, and session minor
    }
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
    bids_subject_id : :obj:`str`
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

def get_file_paths(subject_dir: Path, input_pipeline: str, subject_label: str) -> dict:
    """Get file paths within a subject directory.

    Parameters
    ----------
    subject_dir : :obj:`pathlib.Path`
        The path to the ukb subject directory.
    input_pipeline : :obj:`str`
        The input pipeline used to create the subject directory.
    subject_label : :obj:`str`
        The subject folder identifier.
    

    Returns
    -------
    file_paths : :obj:`dict`
        A dictionary of file paths.
    """

    # Define file path organization for different pipelines
    # Relative to input_dir/subject_dir
    ORGANIZATIONS = {
        "hcpya": {
            "bvals": ["T1w", "Diffusion", "bvals"],
            "bvecs": ["T1w", "Diffusion", "bvecs"],
            "dwi": ["T1w", "Diffusion", "data.nii.gz"],
            "t1w_brain": ["T1w", "T1w_acpc_dc_restore_brain.nii.gz"],
            "brain_mask": ["T1w", "brainmask_fs.nii.gz"],
            "subject2MNI": ["MNINonLinear", "xfms", "acpc_dc2standard.nii.gz"], # Note this is MNI152NLin6Asym
            "freesurfer": ["T1w", subject_label],
        },
        "ukb": {
            "bvals": ["DTI", "dMRI", "dMRI", "bvals"],
            "bvecs": ["DTI", "dMRI", "dMRI", "bvecs"],
            "dwi": ["DTI", "dMRI", "dMRI", "data_ud.nii.gz"],
            "dwiref": ["DTI", "dMRI", "dMRI", "dti_FA.nii.gz"],
            # TODO: Add UKB anatomical paths
            # "t1w_brain": ["T1w", "T1w_acpc_dc_restore_brain.nii.gz"],
            # "brain_mask": ["T1w", "brainmask_fs.nii.gz"],
            # "subject2MNI": ["MNINonLinear", "xfms", "acpc_dc2standard.nii.gz"], # Note this is MNI152NLin6Asym
            # "freesurfer": ["T1w", subject_label], # Note that freesurfer dir is T1w/$SUBJECT_ID, still have to fix this
        }
    }

    # Get organization for input pipeline
    organization = ORGANIZATIONS[input_pipeline]

    # Get and return file paths
    file_paths = {}
    for key, value in organization.items():
        file_paths[key] = subject_dir / Path(*value)
    return file_paths

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