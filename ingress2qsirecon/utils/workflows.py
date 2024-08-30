"""
Nipype Workflows for Ingress2Qsirecon
"""

import os
import shutil
from pathlib import Path

from nipype.pipeline.engine import Workflow

from ingress2qsirecon.utils.interfaces import (
    ConformDwi,
    ConvertWarpfield,
    ExtractB0s,
    NIFTItoH5,
)


def parse_layout(subject_layout):
    # Return the dictionary's values, for dynamic parsing of node output names
    return tuple(subject_layout.values())


def create_single_subject_wf(subject_layout):
    """
    Create a nipype workflow to ingest a single subject.

    This function creates a nipype workflow that takes a dictionary of file paths and
    metadata as input, and outputs a BIDS-formatted directory with the ingested data.

    The workflow consists of the following nodes:

    - ``parse_layout``: a node that takes the input dictionary and extracts the individual
      file paths and metadata.
    - ``conform_dwi``: a node that takes the extracted DWI file paths and metadata, and
      saves them to the BIDS layout.
    - ``create_dwiref``: a node that takes the extracted DWI file paths and metadata, and
      creates a mean b0 image if it does not exist.
    - ``convert_warpfield``: a node that takes the extracted FNIRT warp file paths and
      metadata, and converts them to ITK format.
    - ``nii_to_h5``: a node that takes the converted ITK warp file paths and metadata, and
      saves them as ITK H5 files.

    Parameters
    ----------
    subject_layout : dict
        A dictionary of file paths and metadata for a single subject, from ``create_layout`` function.

    Returns
    -------
    wf : nipype.Workflow
        A nipype workflow that operates on a single subject.
    """
    #### WHY DO I HAVE TO REIMPORT THIS STUFF??
    from nipype import (
        Node,
        Workflow,
    )
    from nipype.interfaces.utility import (
        Function,
        IdentityInterface,
    )

    from ingress2qsirecon.utils.workflows import parse_layout

    ####

    subject_name = subject_layout['subject']

    # Make BIDS subject output folder
    bids_base = subject_layout['bids_base']
    if not os.path.exists(bids_base):
        os.makedirs(Path(bids_base / "anat").resolve())
        os.makedirs(Path(bids_base / "dwi").resolve())

    # Some files, like anatomicals, can just be copied over to the output BIDS directory LPS+ THESE FILES (ConformImage)
    for key in ["t1w_brain", "brain_mask"]:
        bids_key = "bids_" + key
        if key in subject_layout.keys():
            if not os.path.exists(subject_layout[bids_key]):
                shutil.copyfile(subject_layout[key], subject_layout[bids_key])

    # Create single subject workflow
    wf_name = f"ingress2qsirecon_single_subject_{subject_name}_wf"
    wf = Workflow(name=wf_name)

    # Define input node for the single subject workflow
    input_node = Node(IdentityInterface(fields=['subject_layout']), name='input_node')
    input_node.inputs.subject_layout = subject_layout

    # Create node to parse the input dictionary into its individual components
    parse_layout_node = Node(
        Function(
            input_names=['subject_layout'],
            output_names=list(subject_layout.keys()),  # Outputs all fields available in the layout
            function=parse_layout,
        ),
        name='parse_layout_node',
    )

    # Create node to conform DWI and save to BIDS layout
    conform_dwi_node = Node(ConformDwi(), name='conform_dwi')
    wf.connect(
        [
            (input_node, parse_layout_node, [('subject_layout', 'subject_layout')]),
            (
                parse_layout_node,
                conform_dwi_node,
                [
                    ("dwi", "dwi_in_file"),
                    ("bvals", "bval_in_file"),
                    ("bvecs", "bvec_in_file"),
                    ("bids_dwi", "dwi_out_file"),
                    ("bids_bvals", "bval_out_file"),
                    ("bids_bvecs", "bvec_out_file"),
                ],
            ),
        ]
    )
    # Write out .b and .bmtxt (DIPY standard)

    # LPS+ anatomicals

    # If subject does not have DWIREF, run node to extract mean b0
    if "dwiref" not in subject_layout.keys():
        create_dwiref_node = Node(ExtractB0s(), name="create_dwiref")
        wf.connect(
            [
                (
                    parse_layout_node,
                    create_dwiref_node,
                    [("bvals", "bval_file"), ("bids_dwi", "dwi_series"), ("bids_dwiref", "b0_average")],
                )
            ]
        )

    # Convert FNIRT nii warps to ITK nii, then ITK nii to ITK H5, then get to MNI2009cAsym space if needed
    # Start with subject2MNI
    if "subject2MNI" in subject_layout.keys():
        convert_warpfield_node_subject2MNI = Node(ConvertWarpfield(), name="convert_warpfield_subject2MNI")
        convert_warpfield_node_subject2MNI.inputs.itk_out_xfm = str(subject_layout["bids_subject2MNI"]).replace(
            ".h5", ".nii.gz"
        )
        nii_to_h5_node_subject2MNI = Node(NIFTItoH5(), name="nii_to_h5_subject2MNI")
        wf.connect(
            [
                (
                    parse_layout_node,
                    convert_warpfield_node_subject2MNI,
                    [("subject2MNI", "fnirt_in_xfm"), ("MNI_ref", "fnirt_ref_file")],
                ),
                (
                    convert_warpfield_node_subject2MNI,
                    nii_to_h5_node_subject2MNI,
                    [("itk_out_xfm", "xfm_nifti_in")],
                ),
                (
                    parse_layout_node,
                    nii_to_h5_node_subject2MNI,
                    [("bids_subject2MNI", "xfm_h5_out")],
                ),
            ]
        )

    # Then MNI2Subject
    if "MNI2subject" in subject_layout.keys():
        convert_warpfield_node_MNI2subject = Node(ConvertWarpfield(), name="convert_warpfield_MNI2subject")
        convert_warpfield_node_MNI2subject.inputs.itk_out_xfm = str(subject_layout["bids_MNI2subject"]).replace(
            ".h5", ".nii.gz"
        )
        nii_to_h5_node_MNI2subject = Node(NIFTItoH5(), name="nii_to_h5_MNI2subject")
        wf.connect(
            [
                (
                    parse_layout_node,
                    convert_warpfield_node_MNI2subject,
                    [("MNI2subject", "fnirt_in_xfm"), ("MNI_ref", "fnirt_ref_file")],
                ),
                (
                    convert_warpfield_node_MNI2subject,
                    nii_to_h5_node_MNI2subject,
                    [("itk_out_xfm", "xfm_nifti_in")],
                ),
                (
                    parse_layout_node,
                    nii_to_h5_node_MNI2subject,
                    [("bids_MNI2subject", "xfm_h5_out")],
                ),
            ]
        )

    return wf


def create_ingress2qsirecon_wf(layouts, name="ingress2qsirecon_wf", base_dir=os.getcwd()):
    """
    Creates the overall ingress2qsirecon workflow.

    Parameters
    ----------
    layouts : list of dict
        A list of dictionaries, one per subject, from the create_layout function.

    name : str, optional
        The name of the workflow. Default is "ingress2qsirecon_wf".

    base_dir : str, optional
        The base directory in which to create the workflow directory. Default is the current
        working directory.

    Returns
    -------
    wf : nipype.Workflow
        The workflow with the nodes and edges defined.

    """
    wf = Workflow(name=name, base_dir=base_dir)

    subjects_to_run = [layout["subject"] for layout in layouts]
    print(f"Subject(s) to run: {subjects_to_run}")

    for subject_layout in layouts:
        single_subject_wf = create_single_subject_wf(subject_layout)
        wf.add_nodes([single_subject_wf])

    return wf
