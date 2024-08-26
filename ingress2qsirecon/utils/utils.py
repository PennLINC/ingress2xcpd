#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import re
from pathlib import Path
from warnings import warn

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

DIR_PATTERNS = {
    'hcpya': re.compile("(\d+)"),  # HCP is just a number
    'ukb': re.compile("(\d+)_(\d+)_(\d+)"),  # UKB has three numbers for subject, session major, and session minor
}
# List of files that are required for any recon (some pipelines require others too, this is minimum)
required_for_any_recon = [
    "bvals",
    "bvecs",
    "dwi",
]


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
            "subject2MNI": ["MNINonLinear", "xfms", "acpc_dc2standard.nii.gz"],  # Note this is MNI152NLin6Asym
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
        },
    }

    # Get organization for input pipeline
    organization = ORGANIZATIONS[input_pipeline]

    # Get and return file paths
    file_paths = {}
    for key, value in organization.items():
        file_paths[key] = subject_dir / Path(*value)
    return file_paths


class _CreateLayoutInputSpec(BaseInterfaceInputSpec):
    input_dir = Directory(exists=True, mandatory=True, desc="Path to the input directory")
    input_pipeline = Str(mandatory=True, desc="Name of the input pipeline")
    participant_label = traits.List(Str, desc="List of participant labels to search for")


class _CreateLayoutOutputSpec(TraitedSpec):
    layout = traits.List(traits.Dict, desc="List of dictionaries with subject ID, session ID, and paths")


class CreateLayout(SimpleInterface):
    input_spec = _CreateLayoutInputSpec
    output_spec = _CreateLayoutOutputSpec

    def _run_interface(self, runtime):
        input_dir = Path(self.inputs.input_dir)
        input_pipeline = self.inputs.input_pipeline
        participant_label = self.inputs.participant_label
        pattern = DIR_PATTERNS[input_pipeline]  # search pattern for subject ID, session etc

        layout = []
        for potential_dir in input_dir.iterdir():  # loop through all directories in input_dir
            if participant_label and not potential_dir.name.startswith(tuple(participant_label)):
                # Skip if folder in loop does not start with an expected participant label
                continue

            match = re.match(pattern, potential_dir.name)
            if not match:
                # Skip if subject folder does not match expected pattern
                continue

            # Look for data files, make sure at least required ones are present
            # Placeholder for missing_from_subject_dir function
            # if missing_from_subject_dir(potential_dir):
            #     continue

            # If passes all checks, add to layout
            if input_pipeline == "hcpya":
                subject = match.group(1)
                fake_dwi_file = f"/bids/sub-{subject}/dwi/sub-{subject}_dwi.nii.gz"
            elif input_pipeline == "ukb":
                subject, ses_major, ses_minor = match.groups()
                renamed_ses = "%02d%02d" % (int(ses_major), int(ses_minor))
                fake_dwi_file = f"/bids/sub-{subject}/ses-{renamed_ses}/dwi/sub-{subject}_ses-{renamed_ses}_dwi.nii.gz"

            file_paths = get_file_paths(potential_dir, input_pipeline, subject)
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
                "session": None,
                "path": Path(potential_dir),
                "bids_dwi_file": fake_dwi_file,
            }
            # Add file paths to subject layout if they exist
            subject_layout.update({file_type: path for file_type, path in file_paths.items() if os.path.exists(path)})
            layout.append(subject_layout)

        # Sort layout by subject ID
        layout = sorted(layout, key=lambda x: x["subject"])

        # Raise warnings for requested subjects not in layout
        missing_particpants = sorted(set(participant_label) - set([subject["subject"] for subject in layout]))
        if missing_particpants:
            warn(
                f"Requested {missing_particpants} not found in layout, please confirm their data exists and are properly organized."
            )

        # Stop code if layout is empty
        if len(layout) == 0:
            raise ValueError("No subjects found in layout.")

        self._results = {}
        self._results['layout'] = layout
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
