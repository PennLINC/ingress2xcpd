"""

"""
from __future__ import annotations

from beartype import beartype
from .parser import _build_parser
from .main import ingress2qsirecon
@beartype
def main():
    """
    The main entry point of the application.

    This function is responsible for parsing command line arguments and running the main code.

    Parameters:
    None

    Returns:
    None
    """
    # Parse arguments and run the main code
    parser = _build_parser()
    args = parser.parse_args()
    args_dict = vars(args)
    ingress2qsirecon(**args_dict)  

if __name__ == '__main__':
    main()
