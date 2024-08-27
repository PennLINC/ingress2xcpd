from __future__ import annotations

import os
import shutil
from pathlib import Path

from beartype import beartype
from nipype import (
    Node,
    Workflow,
)

from ingress2qsirecon.cli.parser import _build_parser
from ingress2qsirecon.utils.utils import (
    CollectParticipants,
    CreateLayout,
)


@beartype
def ingress2qsirecon(**kwargs):
    # Get the command line arguments
    input_dir = Path(kwargs["input_dir"])
    output_dir = Path(kwargs["output_dir"])
    input_pipeline = kwargs["input_pipeline"]
    participant_label = kwargs["participant_label"]
    work_dir = Path(kwargs["work_dir"])
    check_gradients = kwargs["check_gradients"]
    dry_run = kwargs["dry_run"]
    symlink = kwargs["symlink"]

    # Raise NotImplemented errors for options not implemented yet
    if check_gradients or dry_run or symlink:
        raise NotImplementedError("--check_gradients, --dry_run, and --symlink are not implemented yet.")

    # TMP REMOVE WORK_DIR
    if work_dir.exists():
        shutil.rmtree(work_dir, ignore_errors=True)
    # if workdir doesn't exist, create it
    if not work_dir.exists():
        work_dir.mkdir(parents=True)
    os.chdir(work_dir)

    # If output_dir doesn't exist, create it
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    # Move BIDS scaffold files there
    #ingress2reccon_dir = Path(__file__).parent
    #shutil.copytree(input_dir, output_dir / "scaffold")

    # If participant_label not defined, make it empty list
    if participant_label is None:
        participant_label = []

    # Create overall workflow
    workflow = Workflow(name="ingress2qsirecon_wf", base_dir=work_dir)

    # Create subject gathering nodes
    create_layout_node = Node(CreateLayout(), name="create_layout")
    create_layout_node.inputs.input_dir = input_dir
    create_layout_node.inputs.output_dir = output_dir
    create_layout_node.inputs.input_pipeline = input_pipeline
    create_layout_node.inputs.participant_label = participant_label

    workflow.add_nodes([create_layout_node])
    workflow.run()

    # collect_participants_node = Node(CollectParticipants(), name="collect_participants")
    # collect_participants_node.inputs.participant_label = participant_label

    # workflow.add_nodes([create_layout_node, collect_participants_node])
    # workflow.connect(create_layout_node, "layout", collect_participants_node, "layout")
    # workflow.run()

    # participants_to_run = collect_participants_node.outputs.participants
    # print(f"Participants to run: {participants_to_run}")


def main():
    """
    The main entry point of the CLI.

    This function is responsible for parsing command line arguments and running the main code.
    """
    parser = _build_parser()
    args = parser.parse_args()
    args_dict = vars(args)
    ingress2qsirecon(**args_dict)


if __name__ == '__main__':
    main()
