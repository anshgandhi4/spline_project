import matplotlib.pyplot as plt
import numpy as np

def solve(t, v, endCond, slopes=[]):
    numPts = len(t)
    A = np.ndarray((4 * (numPts - 1), 4 * (numPts - 1)))
    b = np.zeros(4 * (numPts - 1))

    for f in range(numPts - 1):
        stepT = t[f + 1] - t[f]
        A[4 * f] = np.concatenate((np.zeros(4 * f), np.array([0, 0, 0, 1]), np.zeros(4 * (numPts - f - 2))))
        A[4 * f + 1] = np.concatenate((np.zeros(4 * f), np.array([stepT**3, stepT**2, stepT, 1]), np.zeros(4 * (numPts - f - 2))))
        if f <= numPts - 3:
            A[4 * f + 2] = np.concatenate((np.zeros(4 * f), np.array([3 * stepT**2, 2 * stepT, 1, 0, 0, 0, -1, 0]), np.zeros(4 * (numPts - f - 3))))
            A[4 * f + 3] = np.concatenate((np.zeros(4 * f), np.array([6 * stepT, 2, 0, 0, 0, -2, 0, 0]), np.zeros(4 * (numPts - f - 3))))
        b[4 * f] = v[f]
        b[4 * f + 1] = v[f + 1]

    if endCond == "natural first":
        A[-2] = np.concatenate((np.array([0, 0, 1, 0]), np.zeros(4 * (numPts - 2))))
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 2)), np.array([3 * (t[-1] - t[-2])**2, 2 * (t[-1] - t[-2]), 1, 0])))
    elif endCond == "natural second":
        A[-2] = np.concatenate((np.array([0, 2, 0, 0]), np.zeros(4 * (numPts - 2))))
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 2)), np.array([6 * (t[-1] - t[-2]), 2, 0, 0])))
    elif endCond == "periodic":
        A[-2] = np.concatenate((np.array([0, 0, 1, 0]), np.zeros(4 * (numPts - 3)), np.array([-3 * (t[-1] - t[-2])**2, -2 * (t[-1] - t[-2]), -1, 0])))
        A[-1] = np.concatenate((np.array([0, 2, 0, 0]), np.zeros(4 * (numPts - 3)), np.array([-6 * (t[-1] - t[-2]), -2, 0, 0])))
    elif endCond == "not-a-knot":
        A[-2] = np.concatenate((np.array([6, 0, 0, 0, -6, 0, 0, 0]), np.zeros(4 * (numPts - 3))))
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 3)), np.array([6, 0, 0, 0, -6, 0, 0, 0])))
    elif endCond == "specified slope":
        A[-2] = np.concatenate((np.array([0, 0, 1, 0]), np.zeros(4 * (numPts - 2))))
        b[-2] = slopes[0]
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 2)), np.array([3 * (t[-1] - t[-2])**2, 2 * (t[-1] - t[-2]), 1, 0])))
        b[-1] = slopes[1]

    return np.linalg.solve(A, b)

def findInterval(a, t):
    if a == t[-1]:
        return len(t) - 2
    for i in range(len(t)):
        if a >= t[i]:
            continue
        else:
            return i - 1

def f(lst, x, t):
    output = []
    for a in lst:
        i = findInterval(a, t)
        output.append(float(x[4 * i] * (a - t[i])**3 + x[4 * i + 1] * (a - t[i])**2 + x[4 * i + 2] * (a - t[i]) + x[4 * i + 3]))
    return output

def main():
    np.set_printoptions(suppress = True)

    numPts = int(input("Number of points: "))
    t = []
    v = []
    for i in range(numPts):
        t.append(float(input("X" + str(i + 1) + ": ")))
        v.append(float(input("Y" + str(i + 1) + ": ")))

    x = solve(t, v, "natural second")

    a = np.linspace(t[0], t[-1], 1000)
    y = f(a, x, t)

    plt.scatter(t, v, s = 30)
    plt.plot(t, v)
    plt.fill_between(a, y, step = "pre", alpha = 0.4)
    plt.plot(a, y)

    plt.xlabel('Time (seconds)')
    plt.ylabel('Velocity (mph)')
    plt.suptitle('Velocity vs. Time', fontsize = 14)
    # plt.title(f'Distance to Reach {90} MPH: {round(integrate(90, x, t), 4)} miles', fontsize = 10)
    plt.show()

if __name__ == '__main__':
    main()