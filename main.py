import matplotlib.pyplot as plt
import numpy as np

def solve(t, v, extraCond, slopes = []):
    numPts = len(t)
    A = np.ndarray((4 * (numPts - 1), 4 * (numPts - 1)))
    b = np.zeros(4 * (numPts - 1))
    
    for c in range(numPts - 1):
        stepT = t[c + 1] - t[c]
        A[4 * c] = np.concatenate((np.zeros(4 * c), np.array([0, 0, 0, 1]), np.zeros(4 * (numPts - c - 2))))
        A[4 * c + 1] = np.concatenate((np.zeros(4 * c), np.array([stepT**3, stepT**2, stepT, 1]), np.zeros(4 * (numPts - c - 2))))
        if c <= numPts - 3:
            A[4 * c + 2] = np.concatenate((np.zeros(4 * c), np.array([3 * stepT**2, 2 * stepT, 1, 0, 0, 0, -1, 0]), np.zeros(4 * (numPts - c - 3))))
            A[4 * c + 3] = np.concatenate((np.zeros(4 * c), np.array([6 * stepT, 2, 0, 0, 0, -2, 0, 0]), np.zeros(4 * (numPts - c - 3))))
        b[4 * c] = v[c]
        b[4 * c + 1] = v[c + 1]

    if extraCond == "natural first":
        A[-2] = np.concatenate((np.array([0, 0, 1, 0]), np.zeros(4 * (numPts - 2))))
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 2)), np.array([3 * (t[-1] - t[-2])**2, 2 * (t[-1] - t[-2]), 1, 0])))
    elif extraCond == "natural second":
        A[-2] = np.concatenate((np.array([0, 2, 0, 0]), np.zeros(4 * (numPts - 2))))
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 2)), np.array([6 * (t[-1] - t[-2]), 2, 0, 0])))
    elif extraCond == "periodic":
        A[-2] = np.concatenate((np.array([0, 0, 1, 0]), np.zeros(4 * (numPts - 3)), np.array([-3 * (t[-1] - t[-2])**2, -2 * (t[-1] - t[-2]), -1, 0])))
        A[-1] = np.concatenate((np.array([0, 2, 0, 0]), np.zeros(4 * (numPts - 3)), np.array([-6 * (t[-1] - t[-2]), -2, 0, 0])))
    elif extraCond == "not-a-knot":
        A[-2] = np.concatenate((np.array([6, 0, 0, 0, -6, 0, 0, 0]), np.zeros(4 * (numPts - 3))))
        A[-1] = np.concatenate((np.zeros(4 * (numPts - 3)), np.array([6, 0, 0, 0, -6, 0, 0, 0])))
    elif extraCond == "specified slope":
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
        t.append(float(input("T" + str(i + 1) + ": ")))
        v.append(float(input("V" + str(i + 1) + ": ")))

    extraCond = input("Extra condition type (natural first, natural second, periodic, not-a-knot, specified slope): ")
    if extraCond == "specified slope":
        slopes = [input("Start slope: "), input("End slope: ")]
    else:
        slopes = []

    x = solve(t, v, extraCond, slopes)

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