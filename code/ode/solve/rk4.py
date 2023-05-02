from __future__ import annotations

from ..ode import Ode, States
from .solve import Solver, Times


class Rk4(Solver):
    """The Runke-Kutta 4th order solver."""

    def __init__(self, model: Ode, states: States, times: Times) -> None:
        super().__init__(model, states, times)

    def solve_t(self, t: float) -> Rk4:
        constant = 1.0 / 6.0
        one_half = 0.5

        K1 = self._step(self._states, t)
        K2 = self._step(
            states={
                k: s + K1[i] * one_half for i, (k, s) in enumerate(self._states.items())
            },
            t=t + one_half * self._dt,
        )
        K3 = self._step(
            states={
                k: s + K2[i] * one_half for i, (k, s) in enumerate(self._states.items())
            },
            t=t + one_half * self._dt,
        )
        K4 = self._step(
            states={k: s + K3[i] for i, (k, s) in enumerate(self._states.items())},
            t=t + self._dt,
        )

        self.states = {
            k: s + constant * (K1[i] + (K2[i] + K3[i]) * 1 / one_half + K4[i])
            for i, (k, s) in enumerate(self._states.items())
        }

        return self

    def _step(self, states: States, t: float) -> States:
        derivatives = self._model.derivatives(states, t)
        return [self._dt * d for d in derivatives]
