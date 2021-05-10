import numpy as np
import matplotlib.pyplot as plt


def show_pharmacophore(obj, hash=False):
    """
    Demonstration of graph in 3D.

    input(coords, cent):
    'coords' - matrix of coordinates
    'cent' - vector of center mass
    'nodes' - nodes
    """
    if not isinstance(obj, Pharmacophore):
        raise TypeError()
    if not isinstance(hash, bool):
        raise TypeError()
    if hash:
        nodes = [i[0] for i in obj.hash]
        coords = np.array([list(i[1]) for i in obj.hash])
    else:
        nodes = obj.nodes
        coords = obj.coords

    # Assign special color to every type of centers
    colors = {"A": "red", "a": "orange", "D": "plum", "H": "cyan", "P": "olive", "N": "lime"}

    fig = plt.figure(figsize=(30, 24))
    ax = fig.add_subplot(111, projection='3d')

    # Form arrays of xs, ys, and zs from coordinates
    xs = coords[:, 0]
    ys = coords[:, 1]
    zs = coords[:, 2]

    for i in range(len(xs)):
        ax.scatter(xs[i], ys[i], zs[i], c=colors[nodes[i]], marker='o', s=300)

    # Labels on nodes:
    for i in range(len(xs)):
        ax.text(xs[i]+rnd.random()*0.2, ys[i]+rnd.random()*0.2, zs[i]+rnd.random()*0.2, nodes[i], fontsize=20)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    
    # Make simple, bare axis lines through space:
    ax.plot((min(xs) - 1, max(xs) + 1), (0, 0), (0, 0), 'k')
    ax.plot((0, 0), (min(ys) - 1, max(ys) + 1), (0, 0), 'k')
    ax.plot((0, 0), (0, 0), (min(zs) - 1, max(zs) + 1), 'k')

    plt.show()
