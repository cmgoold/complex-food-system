from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Sequence

from ..ode import Ode, States

Times = Sequence[float]

class Solver(ABC):
    """Solve ODE models."""

    def __init__(self, model: Ode, states: States, times: Times) -> None:
        self._model = model
        self._states = states
        self._times = times
        self._dt = times[1] - times[0]
        self._cache: Dict[float, States] = {}

    model = property(lambda self: self._model)
    times = property(lambda self: self._times)
    dt = property(lambda self: self._dt)
    cache = property(lambda self: self._cache)
    
    def solve(self, cache: bool = False) -> Solver:
        for t in self._times:
            if t in self._cache:
                continue
            else:
                self.solve_t(t)
                if cache:
                    self._cache[t] = self._states

        return self

    @abstractmethod
    def solve_t(self, t: float) -> Solver:
        pass

    @property
    def states(self) -> States:
        return self._states

    @states.setter
    def states(self, states) -> None:
        self._states = states


