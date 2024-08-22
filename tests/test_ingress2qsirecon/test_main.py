from __future__ import annotations

from ingress2qsirecon.cli.main import get_hello
from ingress2qsirecon.cli.parser import parse_args


def it_prints_hi_to_the_project_author() -> None:
    expected = 'Hello, Steven Meisler!'
    actual = get_hello('Steven Meisler')
    assert actual == expected
