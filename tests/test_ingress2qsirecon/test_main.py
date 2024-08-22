from __future__ import annotations

from ingress2qsirecon.main import get_hello

def it_prints_hi_to_the_project_author() -> None:
    expected = 'Hello, Steven Meisler!'
    actual = get_hello('Steven Meisler')
    assert actual == expected
