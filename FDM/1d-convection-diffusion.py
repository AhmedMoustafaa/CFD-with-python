import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import imageio

np.set_printoptions(precision=5)


def is_convergent(array1, array2, epsilon):
    diff = np.subtract(np.abs(np.array(array2)), np.abs(np.array(array1)))
    sum1 = np.sum(diff)
    if abs(sum1) < epsilon:
        print(f"error is {sum1}")
        return True
    else:
        print(f"error is {sum1}")
        return False


def conv_diff(n, h, epsilon, boundary, u):
    """solve 1d convection diffusion equations in form :{d^2(f)}/{dx^2} + u*{df}/{dx} using the finite differencing
    method."""
    t1 = np.insert(boundary, 1, np.zeros(n - 2))
    print(t1)
    t_new = t1.copy()
    y = 1
    fig = plt.figure()
    plt.title("Temperature gradient along the x direction")
    writer = PillowWriter(fps=5000)
    frames1 = np.array([t_new])
    plt.title("Temperature gradient along the x direction")
    plt.yticks([])
    while True:
        for i in range(len(t1) - 2):
            t_new[i + 1] = t1[i + 2] * (u * h / 4 + 0.5) + t1[i] * (0.5 - u * h / 4)
        if is_convergent(t1, t_new, epsilon):
            print('will break:')
            print(t_new)
            print('number of iterations: {0}'.format(y))
            break
        else:
            t1 = t_new.copy()
            y = y + 1
            if y%500==0:
                frames2 = np.append(frames1, [t_new], axis=0)
                frames1 = frames2.copy()

    with writer.saving(fig, "animation.gif", 200):
        for frame in frames1:
            data = np.outer(np.ones(10), frame)
            plt.imshow(data, cmap="Reds", alpha=0.8)
            writer.grab_frame()
            plt.pause(0.001)
    return t_new


# example
t2 = conv_diff(200, 1 / 200, 1e-6, np.array([0.0, 1.0]), 0)
