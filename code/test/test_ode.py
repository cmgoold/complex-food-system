import numpy as np

from ode import Logistic, Rk4

def test_logistic_model():
    dt = 0.01
    init = 10
    times = np.arange(0, 1000, dt)
    logistic = Logistic(parameters={"r": 0.1, "k": 100})
    rk4 = Rk4(logistic, states=[init], dt=dt)
    rk4.solve(times)
    assert np.isclose(rk4.states[0], 100, rtol=1e-5)
