import numpy as np

from ode.ode import Logistic
from ode.solve import Rk4


def test_logistic_model():
    dt = 0.0001
    init = 10
    times = np.arange(0, 1+dt, dt)
    logistic = Logistic(parameters={"r": 0.1, "k": 100}, state_names="y")
    rk4 = Rk4(logistic, states={"y": 10})
    rk4.solve(times)
    assert np.isclose(rk4.states["y"], 10.936687, rtol=1e-5)
