import os
import shutil
from pathlib import Path

from nipype.interfaces.utility import (
    Function,
    IdentityInterface,
)
from nipype.pipeline.engine import (
    MapNode,
    Node,
    Workflow,
)

from ingress2qsirecon.utils.interfaces import (
    ConformDwi,
    ConvertWarpfield,
    ExtractB0s,
    NIFTItoH5,
)


def parse_layout(subject_layout):
    # Return the dictionary's values, which will correspond to the dynamic output names of the node
    return tuple(subject_layout.values())


def create_single_subject_wf(subject_layout):

    #### WHY DO I HAVE TO REIMPORT THIS STUFF??
    from nipype import (
        MapNode,
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

    # Some files, like anatomicals, can just be copied over to the output BIDS directory
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
            output_names=list(subject_layout.keys()),  # Set dynamic output names
            function=parse_layout,
        ),
        name='parse_layout_node',
    )

    # Create node to conform DWI and save to BIDS layout
    conform_dwi_node = Node(ConformDwi(), name='conform_dwi')

    # Connect nodes
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

    for subject_layout in layouts:
        single_subject_wf = create_single_subject_wf(subject_layout)
        wf.add_nodes([single_subject_wf])

    # Define input node for the overall workflow
    # input_node = Node(IdentityInterface(fields=['layouts']), name='input_node')

    # Define the MapNode that will run the single subject workflow in parallel
    # single_subject_wf = MapNode(
    #     Function(input_names=['subject_layout'], output_names=['out'], function=create_single_subject_wf),
    #     name='init_ingress2qsirecon_single_subject_wf',
    #     iterfield=['subject_layout'],
    #  )

    # Connect the input node to the parallelized single subject workflows
    # wf.connect([(input_node, single_subject_wf, [('layouts', 'subject_layout')])])

    # Set the layout dictionary list as the input
    # wf.inputs.input_node.layouts = layouts

    return wf
