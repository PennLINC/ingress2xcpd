from __future__ import annotations

from ingress2qsirecon.cli.main import main


def it_prints_hi_to_the_project_author() -> None:
    expected = 'Hello, Steven Meisler!'
    actual = main('Steven Meisler')
    assert actual == expected
