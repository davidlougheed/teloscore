import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["savefig.dpi"] = 300


__all__ = ["plot_versus"]


def plot_versus(file: str):
    data_matrix = []

    s2_hdr = []

    with open(file, "r") as fh:
        s1_hdr = fh.readline().strip().split("\t")

        for line in fh:
            data = line.strip().split("\t")
            s2_hdr.append(data[0])
            data_matrix.append([float(x) for x in data[1:]])

    dm = np.array(data_matrix)

    im = plt.imshow(np.flip(dm, axis=1).T, interpolation="nearest")

    plt.yticks(list(range(len(s1_hdr))), labels=s1_hdr[::-1], fontsize=4)
    plt.xticks(list(range(len(s2_hdr))), labels=s2_hdr, rotation=90, fontsize=4)

    plt.ylabel("Child chromosome index")
    plt.xlabel("Parental chromosome index")

    plt.colorbar(im)

    plt.show()
