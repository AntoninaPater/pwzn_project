import pytest

from project.waveeq_cn import (
    psi_init,
    kin_energy,
    pot_energy,
    v, norm)


def test_kin_energy(arg=psi_init(), expected=0.063109052):
    result = kin_energy(arg)
    assert result == pytest.approx(expected)


def test_pot_energy(arg1=psi_init(), arg2=v(), expected=9.7596373e-10):
    result = pot_energy(arg1, arg2)
    assert result == pytest.approx(expected)


def test_norm(arg=psi_init(), expected=0.9999999):
    result = norm(arg)
    assert result == pytest.approx(expected)