from __future__ import annotations

from ..ode import Ode, States, Parameters, Derivatives, StateNames


class FoodSystem(Ode):
    """The complex food system model."""

    def __init__(self, parameters: Parameters, state_names: StateNames) -> None:
        super().__init__(parameters, state_names)

    a = property(lambda self: self.parameters["a"])
    b = property(lambda self: self.parameters["b"])
    e = property(lambda self: self.parameters["e"])
    f = property(lambda self: self.parameters["f"])
    g = property(lambda self: self.parameters["g"])
    w = property(lambda self: self.parameters["w"])
    s = property(lambda self: self.parameters["s"])
    k = property(lambda self: self.parameters["k"])
    h = property(lambda self: self.parameters["h"])
    m = property(lambda self: self.parameters["m"])
    q = property(lambda self: self.parameters["q"])
    r = property(lambda self: self.parameters["r"])

    def derivatives(self, states: States, t: float) -> Derivatives:
        C, I, D, P = states.values()
        dC_dt = self.a * C * (P / self.b - 1) - self.e * C
        dI_dt = (
            self.f * self.g * C
            - self.w * I
            - D * I / (D * self.s + I)
            + self.k * (self.h - self.f * self.g * C)
        )
        dD_dt = self.m * (self.h * self.q / P - D)
        dP_dt = self.r * P * (self.s * D / I - 1)

        return dC_dt, dI_dt, dD_dt, dP_dt


class DimensionlessFoodSystem(Ode):
    """The dimensionless food system model."""

    def __init__(self, parameters: Parameters, state_names: StateNames) -> None:
        super().__init__(parameters, state_names)

    alpha = property(lambda self: self.parameters["alpha"])
    beta = property(lambda self: self.parameters["beta"])
    delta = property(lambda self: self.parameters["delta"])
    omega = property(lambda self: self.parameters["omega"])
    gamma = property(lambda self: self.parameters["gamma"])
    kappa = property(lambda self: self.parameters["kappa"])
    mu = property(lambda self: self.parameters["mu"])
    rho = property(lambda self: self.parameters["rho"])

    def derivatives(self, states: States, t: float) -> Derivatives:
        v, x, y, z = states.values()
        dv_dtau = v * (self.alpha * z - 1) - self.beta * v
        dx_dtau = (
            self.delta * v
            - self.omega * x
            - self.gamma * x * y / (y + x)
            + self.kappa * (self.gamma - self.delta * v)
        )
        dy_dtau = self.mu * (z**-1 - y)
        dz_dtau = self.rho * z * (y / x - 1)

        return dv_dtau, dx_dtau, dy_dtau, dz_dtau

    @property
    def critical_ratio(self) -> float:
        return (
            self.alpha
            * (self.omega + self.gamma / 2)
            / (self.kappa * self.gamma * (1 + self.beta))
        )

    @property
    def critical_kappa(self) -> float:
        return (
            self.alpha * (self.omega + self.gamma / 2) / (self.gamma * (1 + self.beta))
        )

    @property
    def surplus_ratio(self) -> float:
        return self.kappa * self.critical_ratio

    @classmethod
    def from_dimensional(
        cls, model: FoodSystem, initial_capital: float
    ) -> DimensionlessFoodSystem:
        dimensionless_parameters = {
            "alpha": model.q / model.b,
            "beta": model.e / model.a,
            "delta": model.f
            * model.g
            * initial_capital
            / (model.a * model.h * model.s),
            "omega": model.w / model.a,
            "gamma": (model.a * model.s) ** -1,
            "kappa": model.k,
            "mu": model.m / model.a,
            "rho": model.r / model.a,
        }

        return cls(
            parameters=dimensionless_parameters, state_names=["v", "x", "y", "z"]
        )
