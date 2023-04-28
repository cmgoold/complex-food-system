from abc import ABC, abstractmethod
from typing import Dict, Optional, List, Sequence, Union

States = Dict[str, float]
StateNames = Union[str, List[str]]
Derivatives = Sequence[float]
Parameters = Dict[str, float]

class Ode(ABC):
    """Ordinary differential equation base class."""

    def __init__(self, parameters: Parameters, state_names: StateNames) -> None:
        self._state_names = state_names if isinstance(state_names, list) else [state_names]
        self._n_states = len(self._state_names)
        self._parameters = parameters

    parameters = property(lambda self: self._parameters)
    state_names = property(lambda self: self._state_names)
    n_states = property(lambda self: self._n_states)

    @abstractmethod
    def derivatives(self, states: States, t: float) -> Derivatives:
        if len(states) != self._n_states:
            raise ValueError(f"Number of expected states is {self._n_states} but Ode model was given {len(states)}.")
        if list(states.keys()) == self._state_names:
            raise ValueError(f"`states` dictionary should have key ordering {self._state_names} not {list(states.keys())}.")
        pass
