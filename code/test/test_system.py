import numpy as np

from ode.system import FoodSystem, DimensionlessFoodSystem
from ode.solve import Rk4

BASE_PARAMETERS = {
    "a": 0.05,
    "b": 130,
    "e": 1/(2.5*12),
    "f": 24/12,
    "g": 110*0.75,
    "k": 0.3,
    "h": 130e6,
    "w": 0.3,
    "m": 1/6,
    "q": 160,
    "r": 1/6,
    "s": 1,
}

BASE_STATES = {
    "C": 400e3,
    "I": 130e6,
    "D": 130e6,
    "P": 160,
}



def test_food_system_creation():
    system = FoodSystem(parameters=BASE_PARAMETERS, state_names=list(BASE_STATES.keys()))
    assert isinstance(system, FoodSystem)


def test_dimensionless_food_system_creation():
    system = FoodSystem(parameters=BASE_PARAMETERS, state_names=list(BASE_STATES.keys()))
    dsystem = DimensionlessFoodSystem.from_dimensional(system, 1)
    assert isinstance(dsystem, DimensionlessFoodSystem)
    assert np.isclose(dsystem.critical_ratio, 1.969231, rtol=1e-5)
    assert np.isclose(dsystem.surplus_ratio, 0.5907692, rtol=1e-5)
    assert np.isclose(dsystem.critical_alpha, 0.625, rtol=1e-5)


def test_food_system_results():
    dt = 0.01
    times = np.arange(0, 120, dt)
    system = FoodSystem(parameters=BASE_PARAMETERS, state_names=list(BASE_STATES.keys()))
    rk4 = Rk4(system, states=BASE_STATES, times=times)
    rk4.solve()
    assert np.isclose(list(rk4.states.values()), (328482.7, 96220820, 96129841, 216.2046), rtol=1e-5).all()


def test_dimensionless_food_system_results():
    DSTATES = {
        "v": 1.0, 
        "x": 1.0,
        "y": 1.0,
        "z": 1.0,
    }
    dt_raw = 0.001
    times = np.arange(0, 1 + dt_raw, dt_raw)
    tau = [t / (1 / BASE_PARAMETERS["a"]) for t in times]
    system = FoodSystem(parameters=BASE_PARAMETERS, state_names=list(BASE_STATES.keys()))
    dsystem = DimensionlessFoodSystem.from_dimensional(system, BASE_STATES["C"])
    rk4 = Rk4(dsystem, states=DSTATES, times=tau)
    rk4.solve()
    assert np.isclose(dsystem.critical_ratio, 1.969231, rtol=1e-5)
    assert np.isclose(dsystem.surplus_ratio, 0.5907692, rtol=1e-5)
    assert np.isclose(list(rk4.states.values()), (0.9786494, 0.8860546, 0.9993940, 1.011222), rtol=1e-5).all()
