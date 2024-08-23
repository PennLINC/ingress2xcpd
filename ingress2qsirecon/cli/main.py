from __future__ import annotations

from beartype import beartype

from ingress2qsirecon.cli.parser import _build_parser


@beartype
def ingress2qsirecon(**kwargs):
    print("Hello, World!")
    # Get the command line arguments
    input_dir = kwargs["input_dir"]
    output_dir = kwargs["output_dir"]
    input_pipeline = kwargs["input_pipeline"]
    participant_label = kwargs["participant_label"]
    dry_run = kwargs["dry_run"]
    symlink = kwargs["symlink"]


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
