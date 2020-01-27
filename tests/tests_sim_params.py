import pytest
from requests_mock import mock

from project.sim_params import (
    selector,
    check
    )

test_check_data = [('x0', 41, False),
                   ('l', -2, False),
                   ('sigma', 1, True)]


@pytest.mark.parametrize("arg1, arg2, expected", test_check_data)
def test_check(arg1, arg2, expected):
    result = check(arg1, arg2)
    assert result == expected
