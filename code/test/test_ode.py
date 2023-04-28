import numpy as np

from ode.ode import Logistic, Rk4


def test_logistic_model():
    dt = 0.001
    init = 1
    times = np.arange(0, 1, dt)
    logistic = Logistic(parameters={"r": 0.1, "k": 100}, state_names="y")
    rk4 = Rk4(logistic, states={"y": 10}, dt=dt)
    rk4.solve(times)
    print(rk4.states)
    assert np.isclose(rk4.states["y"], 10.936687, rtol=1e-5)
