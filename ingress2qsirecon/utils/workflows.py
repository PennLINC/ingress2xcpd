"""
Nipype Workflows for Ingress2Qsirecon
"""

import os
from pathlib import Path

from nipype.pipeline.engine import Workflow
from niworkflows.interfaces.images import TemplateDimensions
from templateflow import api as tflow

from ingress2qsirecon.utils.interfaces import (
    Conform,
    ConformDwi,
    ConvertWarpfield,
    ExtractB0s,
    FSLBVecsToTORTOISEBmatrix,
    MRTrixGradientTable,
    NIFTItoH5,
    RobustMNINormalizationRPT,
)


def parse_layout(subject_layout):
    # Return the dictionary's values, for dynamic parsing of node output names
    return tuple(subject_layout.values())


def create_single_subject_wf(subject_layout, skip_mni2009c_norm=False):
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

    # Create single subject workflow
    wf_name = f"ingress2qsirecon_single_subject_{subject_name}_wf"
    wf = Workflow(name=wf_name)

    # Define input node for the single subject workflow
    input_node = Node(
        IdentityInterface(fields=['subject_layout']),
        name='input_node',
    )
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
    # Create node to make b-matrix and bfile from FSL bval/bvec
    create_bmatrix_node = Node(FSLBVecsToTORTOISEBmatrix(), name="create_bmatrix")
    create_bfile_node = Node(MRTrixGradientTable(), name="create_bfile")
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
            (
                conform_dwi_node,
                create_bmatrix_node,
                [
                    ("bval_out_file", "bvals_file"),
                    ("bvec_out_file", "bvecs_file"),
                ],
            ),
            (parse_layout_node, create_bmatrix_node, [("bids_bmtxt", "bmtxt_file")]),
            (
                conform_dwi_node,
                create_bfile_node,
                [
                    ("bval_out_file", "bval_file"),
                    ("bvec_out_file", "bvec_file"),
                ],
            ),
            (parse_layout_node, create_bfile_node, [("bids_b", "b_file_out")]),
        ]
    )

    # Create nodes to conform anatomicals and save to BIDS layout
    # TMP If false because does not work yet
    if "t1w_brain" in subject_layout.keys():
        template_dimensions_node = Node(TemplateDimensions(), name="template_dimensions")
        conform_t1w_node = Node(Conform(), name="conform_t1w")
        wf.connect(
            [
                (
                    parse_layout_node,
                    template_dimensions_node,
                    [("t1w_brain", "t1w_list")],
                ),
                (
                    template_dimensions_node,
                    conform_t1w_node,
                    [("target_shape", "target_shape"), ("target_zooms", "target_zooms")],
                ),
                (
                    parse_layout_node,
                    conform_t1w_node,
                    [("t1w_brain", "in_file"), ("bids_t1w_brain", "out_file")],
                ),
            ]
        )
        if "brain_mask" in subject_layout.keys():
            conform_mask_node = Node(Conform(), name="conform_mask")
            wf.connect(
                [
                    (
                        parse_layout_node,
                        conform_mask_node,
                        [("brain_mask", "in_file"), ("bids_brain_mask", "out_file")],
                    ),
                    (
                        template_dimensions_node,
                        conform_mask_node,
                        [("target_shape", "target_shape"), ("target_zooms", "target_zooms")],
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

    # Convert FNIRT nii warps to ITK nii, then ITK nii to ITK H5
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

    # Now get transform to MNI2009cAsym
    MNI_template = subject_layout["MNI_template"]
    if MNI_template != "MNI152NLin2009cAsym" and skip_mni2009c_norm == False:
        # Get MNI brain and mask
        MNI2009cAsym_brain_path = str(
            tflow.get('MNI152NLin2009cAsym', desc="brain", suffix="T1w", resolution=1, extension=".nii.gz")
        )
        MNI2009cAsym_mask_path = str(
            tflow.get('MNI152NLin2009cAsym', desc="brain", suffix="mask", resolution=1, extension=".nii.gz")
        )
        # Create transform node
        anat_norm_interface = RobustMNINormalizationRPT(float=True, generate_report=True, flavor="precise")
        anat_nlin_normalization = Node(anat_norm_interface, name="anat_nlin_normalization")
        # Set inputs
        anat_nlin_normalization.inputs.template = MNI2009cAsym_brain_path
        anat_nlin_normalization.inputs.reference_image = MNI2009cAsym_brain_path
        anat_nlin_normalization.inputs.reference_mask = MNI2009cAsym_mask_path
        anat_nlin_normalization.inputs.orientation = "LPS"

        # Create output node to save out relevant files
        def save_xfm_outputs(
            to_template_nonlinear_transform_in,
            from_template_nonlinear_transform_in,
            to_template_nonlinear_transform_out,
            from_template_nonlinear_transform_out,
        ):
            import shutil

            # Dictionary of inputs to save
            files = {
                to_template_nonlinear_transform_out: to_template_nonlinear_transform_in,
                from_template_nonlinear_transform_out: from_template_nonlinear_transform_in,
            }

            for filename, file_content in files.items():
                # Copy or move files to the output directory
                shutil.copy(file_content, filename)

        save_outputs_node = Node(
            Function(
                input_names=[
                    "to_template_nonlinear_transform_in",
                    "from_template_nonlinear_transform_in",
                    "to_template_nonlinear_transform_out",
                    "from_template_nonlinear_transform_out",
                ],
                function=save_xfm_outputs,
            ),
            name="save_outputs_node",
        )
        save_outputs_node.inputs.to_template_nonlinear_transform_out = str(subject_layout["bids_subject2MNI"]).replace(
            subject_layout["MNI_template"], "MNI152NLin2009cAsym"
        )
        save_outputs_node.inputs.from_template_nonlinear_transform_out = str(
            subject_layout["bids_MNI2subject"]
        ).replace(subject_layout["MNI_template"], "MNI152NLin2009cAsym")
        # Link T1w brain and mask to node
        wf.connect(
            [
                (
                    conform_t1w_node,
                    anat_nlin_normalization,
                    [("out_file", "moving_image")],
                ),
                (
                    conform_mask_node,
                    anat_nlin_normalization,
                    [("out_file", "moving_mask")],
                ),
                (
                    anat_nlin_normalization,
                    save_outputs_node,
                    [
                        ('composite_transform', 'to_template_nonlinear_transform_in'),
                        ('inverse_composite_transform', 'from_template_nonlinear_transform_in'),
                    ],
                ),
            ]
        )

    return wf


def create_ingress2qsirecon_wf(layouts, name="ingress2qsirecon_wf", base_dir=os.getcwd(), skip_mni2009c_norm=False):
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
        single_subject_wf = create_single_subject_wf(subject_layout, skip_mni2009c_norm=skip_mni2009c_norm)
        wf.add_nodes([single_subject_wf])

    return wf
