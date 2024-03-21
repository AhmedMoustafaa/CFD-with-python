import numpy as np
import matplotlib.pyplot as plt


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


def conv_diff(n,h, epsilon, boundary, u):
    """solve 1d convection diffusion equations in form :{d^2(f)}/{dx^2} + u*{df}/{dx} using the finite differencing
    method."""
    t1 = np.insert(boundary, 1, np.zeros(n-2))
    print(t1)
    t2 = t1.copy()
    y = 1
    while True:
        for i in range(len(t1) - 2):
            t2[i + 1] = t1[i + 2] * (u*h/4 + 0.5) + t1[i] * (0.5 - u*h/4)
        if is_convergent(t1, t2, epsilon):
            print('will break:')
            print(t2)
            print('number of iterations: {0}'.format(y))
            break
        else:
            t1 = t2.copy()
            y = y + 1
    x = np.arange(n)*h
    plt.plot(x, t2, linewidth=2)
    plt.xlabel('x')
    plt.ylabel('Temperature')
    plt.show()
    # Sample data (replace with your actual data)
    data = np.outer(np.ones(10), t2)  # 10x10 array of random values between 0 and 1

    # Create a heatmap with a red colormap
    plt.imshow(data, cmap="Reds", alpha=0.9)

    # Add labels and title
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Data Visualization with Reds Colormap")
    plt.show()
    return t2

# example
t2 = conv_diff(100, 1/100, 1e-6, np.array([0.0, 1.0]), 2)
