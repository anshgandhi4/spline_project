import matplotlib.pyplot as plt
import numpy as np

def solve(t, v):
    A = np.matrix([[0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [(t[1] - t[0])**3, (t[1] - t[0])**2, (t[1] - t[0]), 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, (t[2] - t[1])**3, (t[2] - t[1])**2, (t[2] - t[1]), 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0, 0, 0, (t[3] - t[2])**3, (t[3] - t[2])**2, (t[3] - t[2]), 1],
                    [3*(t[1] - t[0])**2, 2*(t[1] - t[0]), 1, 0, 0, 0, -1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 3*(t[2] - t[1])**2, 2*(t[2] - t[1]), 1, 0, 0, 0, -1, 0],
                    [3*(t[1] - t[0]), 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 3*(t[2] - t[1]), 1, 0, 0, 0, -1, 0, 0],
                    [0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 3*(t[3] - t[2]), 1, 0, 0]])

    b = np.array([[v[0]], [v[1]], [v[1]], [v[2]], [v[2]], [v[3]], [0], [0], [0], [0], [0], [0]])
    
    return np.linalg.solve(A, b)

def f(lst, x, t):
    output = []
    for a in lst:
        if t[0] <= a < t[1]:
            output.append(float(x[:4][0] * (a - t[0])**3 + x[:4][1] * (a - t[0])**2 + x[:4][2] * (a - t[0]) + x[:4][3]))
        elif t[1] <= a < t[2]:
            output.append(float(x[4:8][0] * (a - t[1])**3 + x[4:8][1] * (a - t[1])**2 + x[4:8][2] * (a - t[1]) + x[4:8][3]))
        elif t[2] <= a <= t[3]:
            output.append(float(x[8:][0] * (a - t[2])**3 + x[8:][1] * (a - t[2])**2 + x[8:][2] * (a - t[2]) + x[8:][3]))
        else:
            output.append(-1.0)
    return output

def main():
    np.set_printoptions(suppress = True)

    t = [0, 3.1, 10.3, 30.1]
    v = [0, 30, 60, 90]
    x = solve(t, v)

    a = np.linspace(t[0], t[3], 1000)
    y = f(a, x, t)

    plt.plot(t, v)
    plt.plot(a, y)
    plt.scatter(t, v, s = 30)

    plt.xlabel('Time (seconds)')
    plt.ylabel('Velocity (mph)')
    plt.show()

if __name__ == '__main__':
    main()
