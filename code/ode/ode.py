from abc import ABC, abstractmethod
from typing import Dict, Optional, List, Sequence

States = Sequence[float]
Derivatives = Sequence[float]
Parameters = Dict[str, float]

class Ode(ABC):
    """Ordinary differential equation base class."""

    def __init__(self, parameters: Parameters, n_states: Optional[int] = None) -> None:
        self._n_states = 1 if n_states is None else n_states
        self._parameters = parameters

    parameters = property(lambda self: self._parameters)
    n_states = property(lambda self: self._n_states)

    @abstractmethod
    def derivatives(self, states: States, t: float) -> Derivatives:
        if len(states) != self._n_states:
            raise ValueError(f"Number of expected states is {self._n_states} but Ode model was given {len(states)}.")
        pass

class Logistic(Ode):

    def __init__(self, parameters: Parameters) -> None:
        super().__init__(parameters)

    r = property(lambda self: self._parameters["r"])
    k = property(lambda self: self._parameters["k"])

    def derivatives(self, states: States, t: float) -> Derivatives:
        y, = states
        return self.r * y * (1 - y / self.k),
