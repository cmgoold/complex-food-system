from .ode import Ode, Parameters, States, Derivatives, StateNames

class Logistic(Ode):

    def __init__(self, parameters: Parameters, state_names: StateNames) -> None:
        super().__init__(parameters, state_names)

    r = property(lambda self: self._parameters["r"])
    k = property(lambda self: self._parameters["k"])

    def derivatives(self, states: States, t: float) -> Derivatives:
        y, = states.values()
        return self.r * y * (1 - y / self.k),
