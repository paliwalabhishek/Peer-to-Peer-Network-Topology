import math
import random
from random import randrange
import matplotlib.pyplot as plt
import sys


class node:
    def __init__(self, nid, theta, x, y):
        self.nid = nid
        self.theta = theta
        self.neighbors = []
        self.x = x
        self.y = y


class ring_topology:
    def __init__(self, N, k):
        nodes = []
        for i in range(N):
            theta = ((2 * math.pi) / N) * i
            x = math.cos(theta)
            y = math.sin(theta)
            n = node(i, theta, x, y)
            nodes.append(n)
        for i in range(N):
            ran = random.sample(nodes, len(nodes))
            for j in range(k):
                nodes[i].neighbors.append(ran[j])
        self.compute(N, k, nodes)

    def assigDist(self, N, n1, n2):
        return min(N - abs(n1 - n2), abs(n1 - n2))

    def sumOfDistance(self, N, k, nodes):
        neigDistance = 0
        for i in range(N):
            for j in range(k):
                neigDistance += self.assigDist(N, nodes[i].nid, nodes[i].neighbors[j].nid)
        return neigDistance

    def compute(self, N, k, nodes):
        dist_list = []
        for z in range(40):
            for j in range(N):
                n1 = nodes[j].nid
                index = randrange(k)
                n2 = nodes[j].neighbors[index].nid
                slist = []
                rlist = []
                for m in range(k):
                    slist.append(nodes[n2].neighbors[m].nid)
                    rlist.append(nodes[n2].neighbors[m].nid)
                for m in range(k):
                    slist.append(nodes[n1].neighbors[m].nid)
                    rlist.append(nodes[n1].neighbors[m].nid)
                distSlist = [self.assigDist(N, n1, neigh) for neigh in slist]
                distRlist = [self.assigDist(N, n2, neigh) for neigh in rlist]
                for f in range((2 * k) - 1):
                    for b in range((2 * k) - 1 - f):
                        if distSlist[b] > distSlist[b + 1]:
                            swap = distSlist[b]
                            distSlist[b] = distSlist[b + 1]
                            distSlist[b + 1] = swap
                            swap = slist[b]
                            slist[b] = slist[b + 1]
                            slist[b + 1] = swap
                        if distRlist[b] > distRlist[b + 1]:
                            swap = distRlist[b]
                            distRlist[b] = distRlist[b + 1]
                            distRlist[b + 1] = swap
                            swap = rlist[b]
                            rlist[b] = rlist[b + 1]
                            rlist[b + 1] = swap
                a = 0
                i = 0
                sflag = 0
                while a < k and i < (2 * k) - 1:
                    if a == 0 and slist[i] != n1:
                        nodes[n1].neighbors[a] = nodes[slist[i]]
                        i += 1
                        a += 1
                    for p in range(a):
                        if slist[i] == nodes[n1].neighbors[p] and a != 1:
                            sflag = 1
                    if sflag == 1 and i < (2 * k):
                        i += 1
                    if slist[i] != n1 and slist[i] != nodes[n1].neighbors[a - 1].nid:
                        nodes[n1].neighbors[a] = nodes[slist[i]]
                        a += 1
                        i += 1
                    else:
                        i += 1
                a = 0
                i = 0
                rflag = 0
                while a < k and i < (2 * k) - 1:
                    if a == 0 and rlist[i] != n2:
                        nodes[n2].neighbors[a] = nodes[rlist[i]]
                    for p in range(a):
                        if rlist[i] == nodes[n2].neighbors[p].nid and a != 1:
                            rflag = 1
                    if rflag == 1 and i < (2 * k):
                        i += 1
                    if rlist[i] != n2 and rlist[i] != nodes[n2].neighbors[a - 1].nid:
                        nodes[n2].neighbors[a] = nodes[rlist[i]]
                        a += 1
                        i += 1
                    else:
                        i += 1
            dist = self.sumOfDistance(N, k, nodes)
            dist_list.append(dist)
            if z == 1 or (z % 5 == 0 and z <= 15):
                self.plotTopology(N, k, nodes, z)
                self.printNe(N, k, nodes, z)

            # print(f"Sum Of Distance after {z} cycle : {dist}")
        file = open(f"R_N{N}_k{k}.txt", "w")
        for e in dist_list:
            file.write(str(e) + "\n")
        file.close()
        self.plotDistance(40, dist_list)

    def printNe(self, N, k, nodes, cycle):
        file = open(f"R_N{N}_k{k}_{cycle}.txt", "w")
        for i in range(N):
            file.write(f"Neighbor List of {i}th node is: \n")
            for j in range(len(nodes[i].neighbors)):
                file.write(f"{nodes[i].neighbors[j].nid} \t")
            file.write("\n")
        file.close()

    def plotTopology(self, N, k, nodes, cycle):
        for i in range(N):
            xi = nodes[i].x
            yi = nodes[i].y
            for j in range(2):
                xni = nodes[i].neighbors[j].x
                yni = nodes[i].neighbors[j].y
                plt.plot([xi, xni], [yi, yni], color='r')
        plt.title(f"{cycle} Cycle")
        plt.savefig(f"R_N{N}_k{k}_{cycle}.png")
        plt.show()

    def plotDistance(self, cycle, dist_list):
        plt.plot(range(cycle), dist_list)
        plt.ylabel("Distance")
        plt.xlabel("Cycle")
        plt.show()


class dynamic_ring_topology:
    def __init__(self, N, k, n, r):
        nodes = []
        for i in range(N):
            theta = ((2 * math.pi) / N) * i
            x = math.cos(theta)
            y = math.sin(theta)
            no = node(i, theta, x, y)
            nodes.append(no)
        for i in range(N):
            ran = random.sample(nodes, len(nodes))
            for j in range(k):
                nodes[i].neighbors.append(ran[j])
        self.compute(N, k, n, r, nodes)

    def distance(self, n1, n2, r, nodes):
        x1 = r * (nodes[n1].x)
        y1 = r * (nodes[n1].y)
        x2 = r * (nodes[n2].x)
        y2 = r * (nodes[n2].y)
        return math.sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))

    def sumOfDistance(self, N, k, r, nodes):
        neigDistance = 0
        for i in range(N):
            for j in range(k):
                neigDistance += self.distance(nodes[i].nid, nodes[i].neighbors[j].nid, r, nodes)
        return neigDistance

    def compute(self, N, k, n, r, nodes):
        dist_list = []
        ring_count = 0
        des_r = r[ring_count + 1]
        radius = des_r
        for z in range(40):
            if z % 3 == 0 and radius < des_r:
                radius += 1
            if z % 5 == 0 and ring_count < n:
                des_r = r[ring_count + 1]
            for j in range(N):
                n1 = nodes[j].nid
                index = randrange(k)
                n2 = nodes[j].neighbors[index].nid
                slist = []
                rlist = []
                for m in range(k):
                    slist.append(nodes[n2].neighbors[m].nid)
                    rlist.append(nodes[n2].neighbors[m].nid)
                for m in range(k):
                    slist.append(nodes[n1].neighbors[m].nid)
                    rlist.append(nodes[n1].neighbors[m].nid)
                distSlist = [self.distance(n1, neigh, radius, nodes) for neigh in slist]
                distRlist = [self.distance(n2, neigh, radius, nodes) for neigh in rlist]
                for f in range((2 * k) - 1):
                    for b in range((2 * k) - 1 - f):
                        if distSlist[b] > distSlist[b + 1]:
                            swap = distSlist[b]
                            distSlist[b] = distSlist[b + 1]
                            distSlist[b + 1] = swap
                            swap = slist[b]
                            slist[b] = slist[b + 1]
                            slist[b + 1] = swap
                        if distRlist[b] > distRlist[b + 1]:
                            swap = distRlist[b]
                            distRlist[b] = distRlist[b + 1]
                            distRlist[b + 1] = swap
                            swap = rlist[b]
                            rlist[b] = rlist[b + 1]
                            rlist[b + 1] = swap
                a = 0
                i = 0
                sflag = 0
                while a < k and i < (2 * k) - 1:
                    if a == 0 and slist[i] != n1:
                        nodes[n1].neighbors[a] = nodes[slist[i]]
                        i += 1
                        a += 1
                    for p in range(a):
                        if slist[i] == nodes[n1].neighbors[p] and a != 1:
                            sflag = 1
                    if sflag == 1 and i < (2 * k):
                        i += 1
                    if slist[i] != n1 and slist[i] != nodes[n1].neighbors[a - 1].nid:
                        nodes[n1].neighbors[a] = nodes[slist[i]]
                        a += 1
                        i += 1
                    else:
                        i += 1
                a = 0
                i = 0
                rflag = 0
                while a < k and i < (2 * k) - 1:
                    if a == 0 and rlist[i] != n2:
                        nodes[n2].neighbors[a] = nodes[rlist[i]]
                    for p in range(a):
                        if rlist[i] == nodes[n2].neighbors[p].nid and a != 1:
                            rflag = 1
                    if rflag == 1 and i < (2 * k):
                        i += 1
                    if rlist[i] != n2 and rlist[i] != nodes[n2].neighbors[a - 1].nid:
                        nodes[n2].neighbors[a] = nodes[rlist[i]]
                        a += 1
                        i += 1
                    else:
                        i += 1
            dist = self.sumOfDistance(N, k, radius, nodes)
            dist_list.append(dist)
            if z == 1 or (z % 5 == 0 and z <= 15):
                self.plotTopology(N, k, nodes, z)
                self.printNe(N, k, nodes, z)

        file = open(f"D_N{N}_k{k}.txt", "w")
        for e in dist_list:
            file.write(str(e) + "\n")
        file.close()
        self.plotDistance(40, dist_list)

    def printNe(self, N, k, nodes, cycle):
        file = open(f"D_N{N}_k{k}_{cycle}.txt", "w")
        for i in range(N):
            file.write(f"Neighbor List of {i}th node is: \n")
            for j in range(len(nodes[i].neighbors)):
                file.write(f"{nodes[i].neighbors[j].nid} \t")
            file.write("\n")
        file.close()

    def plotTopology(self, N, k, nodes, cycle):
        for i in range(N):
            xi = nodes[i].x
            yi = nodes[i].y
            for j in range(2):
                xni = nodes[i].neighbors[j].x
                yni = nodes[i].neighbors[j].y
                plt.plot([xi, xni], [yi, yni], color='r')
        plt.title(f"{cycle} Cycle")
        plt.savefig(f"D_N{N}_k{k}_{cycle}.png")
        plt.show()

    def plotDistance(self, cycle, dist_list):
        plt.plot(range(cycle), dist_list)
        plt.ylabel("Distance")
        plt.xlabel("Cycle")
        plt.show()


class cov_topology:
    def __init__(self, N, k):
        nodes = []
        count = -1
        theta_count = 0
        diff_theta = 22.5
        spike_len = 1
        flag = 0
        for i in range(N):
            if count < i < count + 75 and flag == 0:
                theta = math.radians(theta_count * diff_theta)
                x = math.cos(theta) * spike_len
                y = math.sin(theta) * spike_len
                spike_len += 0.01
                n = node(i, theta, x, y)
                nodes.append(n)
                if i + 1 == count + 75:
                    theta_count += 1
                    spike_len = 1
                    count += 74
                    flag = 1
            if count < i < count + 45 and flag == 1:
                theta = math.radians(theta_count * diff_theta)
                x = math.cos(theta) * spike_len
                y = math.sin(theta) * spike_len
                spike_len += 0.01
                n = node(i, theta, x, y)
                nodes.append(n)
                if i + 1 == count + 45:
                    theta_count += 1
                    spike_len = 1
                    count += 44
                    flag = 0
        for i in range(N):
            ran = random.sample(nodes, len(nodes))
            for j in range(k):
                nodes[i].neighbors.append(ran[j])
        self.compute(N, k, nodes)

    def distance(self, n1, n2, nodes):
        x1 = nodes[n1].x
        y1 = nodes[n1].y
        x2 = nodes[n2].x
        y2 = nodes[n2].y
        return math.sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))

    def sumOfDistance(self, N, k, nodes):
        neigDistance = 0
        for i in range(N):
            for j in range(k):
                neigDistance += self.distance(nodes[i].nid, nodes[i].neighbors[j].nid, nodes)
        return neigDistance

    def compute(self, N, k, nodes):
        dist_list = []
        for z in range(40):
            for j in range(N):
                n1 = nodes[j].nid
                index = randrange(k)
                n2 = nodes[j].neighbors[index].nid
                slist = []
                rlist = []
                for m in range(k):
                    slist.append(nodes[n2].neighbors[m].nid)
                    rlist.append(nodes[n2].neighbors[m].nid)
                for m in range(k):
                    slist.append(nodes[n1].neighbors[m].nid)
                    rlist.append(nodes[n1].neighbors[m].nid)
                distSlist = [self.distance(n1, neigh, nodes) for neigh in slist]
                distRlist = [self.distance(n2, neigh, nodes) for neigh in rlist]
                for f in range((2 * k) - 1):
                    for b in range((2 * k) - 1 - f):
                        if distSlist[b] > distSlist[b + 1]:
                            swap = distSlist[b]
                            distSlist[b] = distSlist[b + 1]
                            distSlist[b + 1] = swap
                            swap = slist[b]
                            slist[b] = slist[b + 1]
                            slist[b + 1] = swap
                        if distRlist[b] > distRlist[b + 1]:
                            swap = distRlist[b]
                            distRlist[b] = distRlist[b + 1]
                            distRlist[b + 1] = swap
                            swap = rlist[b]
                            rlist[b] = rlist[b + 1]
                            rlist[b + 1] = swap
                a = 0
                i = 0
                sflag = 0
                while a < k and i < (2 * k) - 1:
                    if a == 0 and slist[i] != n1:
                        nodes[n1].neighbors[a] = nodes[slist[i]]
                        i += 1
                        a += 1
                    for p in range(a):
                        if slist[i] == nodes[n1].neighbors[p] and a != 1:
                            sflag = 1
                    if sflag == 1 and i < (2 * k):
                        i += 1
                    if slist[i] != n1 and slist[i] != nodes[n1].neighbors[a - 1].nid:
                        nodes[n1].neighbors[a] = nodes[slist[i]]
                        a += 1
                        i += 1
                    else:
                        i += 1
                a = 0
                i = 0
                rflag = 0
                while a < k and i < (2 * k) - 1:
                    if a == 0 and rlist[i] != n2:
                        nodes[n2].neighbors[a] = nodes[rlist[i]]
                    for p in range(a):
                        if rlist[i] == nodes[n2].neighbors[p].nid and a != 1:
                            rflag = 1
                    if rflag == 1 and i < (2 * k):
                        i += 1
                    if rlist[i] != n2 and rlist[i] != nodes[n2].neighbors[a - 1].nid:
                        nodes[n2].neighbors[a] = nodes[rlist[i]]
                        a += 1
                        i += 1
                    else:
                        i += 1
            dist = self.sumOfDistance(N, k, nodes)
            dist_list.append(dist)
            if z == 1 or (z % 5 == 0 and z <= 15):
                self.plotTopology(N, k, nodes, z)
                self.printNe(N, k, nodes, z)
        file = open(f"C_N{N}_k{k}.txt", "w")
        for e in dist_list:
            file.write(str(e) + "\n")
        file.close()
        self.plotDistance(40, dist_list)

    def printNe(self, N, k, nodes, cycle):
        file = open(f"C_N{N}_k{k}_{cycle}.txt", "w")
        for i in range(N):
            file.write(f"Neighbor List of {i}th node is: \n")
            for j in range(len(nodes[i].neighbors)):
                file.write(f"{nodes[i].neighbors[j].nid} \t")
            file.write("\n")
        file.close()

    def plotDistance(self, cycle, dist_list):
        plt.plot(range(cycle), dist_list)
        plt.ylabel("Distance")
        plt.xlabel("Cycle")
        plt.show()

    def plotTopology(self, N, k, nodes, cycle):
        x = []
        y = []
        for i in range(17):
            theta = ((2 * math.pi) / 16) * i
            x.append(math.cos(theta))
            y.append(math.sin(theta))
        plt.plot(x, y, color='r')
        for i in range(N):
            xi = nodes[i].x
            yi = nodes[i].y
            for j in range(2):
                xni = nodes[i].neighbors[j].x
                yni = nodes[i].neighbors[j].y
                plt.plot([xi, xni], [yi, yni], color='r')
        plt.title(f"{cycle} Cycle for K value: {k} ")
        plt.savefig(f"C_N{N}_k{k}_{cycle}.png")
        plt.show()


def main(name_of_file, nodes, neig, topology, n, r):
    if topology == 'R':
        ring_topology(nodes, neig)
    if topology == 'D':
        dynamic_ring_topology(nodes, neig, n, r)
    if topology == 'C':
        cov_topology(nodes, neig)


if __name__ == "__main__":
    name_of_file = sys.argv[0]
    nodes = int(sys.argv[1])
    neig = int(sys.argv[2])
    topology = sys.argv[3]
    count_of_radius = 0
    radius = []
    if topology == 'D':
    	count_of_radius = int(sys.argv[4])
    	radius = list(map(int, sys.argv[5:]))
    main(name_of_file, nodes, neig, topology, count_of_radius, radius)
"""
N = 1000
k = 40
n = 5
r = [4, 6, 8, 10, 12]
# ring_topology(N, k)
# dynamic_ring_topology(N, k, n, r)
cov_topology(N, k)
"""