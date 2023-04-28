from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from ..ode import Ode, States


class Solver(ABC):
    """Solve ODE models."""

    def __init__(self, model: Ode, states: States) -> None:
        self._model = model
        self._states = states
        self._dt = None

    model = property(lambda self: self._model)
    
    @property
    def dt(self) -> float:
        return self._dt

    @dt.setter
    def dt(self) -> None:
        self._dt = dt

    def solve(self, times: Sequence[float]) -> Solver:
        self._dt = times[1] - times[0]
        for t in times:
            self.solve_t(t)

        return self

    @abstractmethod
    def solve_t(self, t: float) -> Solver:
        if self._dt is None:
            raise ValueError("dt interval for the solver is `None` -- set the dt interval first.")
        pass

    @property
    def states(self) -> States:
        return self._states

    @states.setter
    def states(self, states) -> None:
        self._states = states


