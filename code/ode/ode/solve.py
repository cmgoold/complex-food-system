from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from .ode import Ode, States


class Solver(ABC):
    """Solve ODE models."""

    def __init__(self, model: Ode, states: States, dt: float) -> None:
        self._model = model
        self._states = states
        self._dt = dt

    model = property(lambda self: self._model)
    dt = property(lambda self: self._dt)

    def solve(self, times: Sequence[float]) -> None:
        for t in times:
            self.solve_t(t)

    @abstractmethod
    def solve_t(self, t: float) -> Solver:
        pass

    @property
    def states(self) -> States:
        return self._states

    @states.setter
    def states(self, states) -> None:
        self._states = states


class Rk4(Solver):
    """The Runke-Kutta 4th order solver."""

    def __init__(self, model: Ode, states: States, dt: float) -> None:
        super().__init__(model, states, dt)

    def solve_t(self, t: float) -> Rk4:
        constant = 1.0 / 6.0
        one_half = 0.5

        K1 = self._euler(self._states, t)
        K2 = self._euler({k: s + K1[i] * one_half for i, (k, s) in enumerate(self._states.items())}, t + one_half * self._dt)
        K3 = self._euler({k: s + K2[i] * one_half for i, (k, s) in enumerate(self._states.items())}, t + one_half * self._dt)
        K4 = self._euler({k: s + K3[i] for i, (k, s) in enumerate(self._states.items())}, t + self._dt)

        self.states = {
            k: s + constant * (K1[i] + (K2[i] + K3[i]) * 1/one_half + K4[i]) 
            for i, (k, s) 
            in enumerate(self._states.items())
        }

        return self

    def _euler(self, states: States, t: float) -> States:
        derivatives = self._model.derivatives(states, t)
        return [self._dt * d for d in derivatives]
