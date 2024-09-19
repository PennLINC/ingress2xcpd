from __future__ import annotations

import os
import shutil
from pathlib import Path

from beartype import beartype

from ingress2qsirecon.cli.parser import _build_parser
from ingress2qsirecon.data import load_resource
from ingress2qsirecon.utils.functions import create_layout
from ingress2qsirecon.utils.workflows import create_ingress2qsirecon_wf


@beartype
def _ingress2qsirecon(**kwargs):
    """
    The main function

    This initializes the import directories, then creates and run the nipype workflow.
    """
    # Get the command line arguments
    input_dir = Path(kwargs["input_dir"])
    output_dir = Path(kwargs["output_dir"])
    input_pipeline = kwargs["input_pipeline"]
    participant_label = kwargs["participant_label"]
    work_dir = Path(kwargs["work_dir"])
    skip_mni2009c_norm = kwargs["skip_mni2009c_norm"]
    check_gradients = kwargs["check_gradients"]
    dry_run = kwargs["dry_run"]
    symlink = kwargs["symlink"]

    # Raise NotImplemented errors for options not implemented yet
    if check_gradients or dry_run or symlink:
        raise NotImplementedError("--check_gradients, --dry_run, and --symlink are not implemented yet.")

    # Create working directory
    if not work_dir.exists():
        work_dir.mkdir(parents=True)
    os.chdir(work_dir)

    # If output_dir doesn't exist, create it
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    # Move BIDS scaffold files to output directory
    bids_scaffold = str(load_resource('bids_scaffold/'))
    if not os.path.exists(os.path.join(output_dir, "dataset_description.json")):
        shutil.copytree(bids_scaffold, output_dir, dirs_exist_ok=True)

    # If participant_label not defined, make it empty list
    if participant_label is None:
        participant_label = []

    # Make list of dictionaries with all information about the ingressions, one dict per subject
    layouts = create_layout(input_dir, output_dir, input_pipeline, participant_label)

    # Create and run overall workflow, which will be broken down to single subject workflows
    ingress2qsirecon_wf = create_ingress2qsirecon_wf(layouts, base_dir=work_dir, skip_mni2009c_norm=skip_mni2009c_norm)
    ingress2qsirecon_wf.run()


def main():
    """
    The main entry point of the CLI.

    This function is responsible for parsing command line arguments and running the main code.
    """
    parser = _build_parser()
    args = parser.parse_args()
    args_dict = vars(args)
    _ingress2qsirecon(**args_dict)


if __name__ == '__main__':
    main()
