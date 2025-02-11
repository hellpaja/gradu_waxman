import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random
import math

class Waxman_square:
    def __init__(self, n, R, gamma = 0.2, beta = 1.0):
        # Alfa and beta are constants on the interval from 0 to 1.
        # Constant n is the number of nodes, R is the length of the side of the square (in kilometers)
        # gamma is the fibers loss rate?. Beta is usually 1 and state-of-the-art gamma is 0.2 dB per km
        # for a state of the art capel alfa*sqrt(2)*R = 226 (km).
        self.beta = beta
        self.n = n
        self.R = R
        self.gamma = gamma
        self.alfa = 226/(np.sqrt(2)*R)

        # alustetaan tarvittavat matriisit nolla matriiseiksi, joissa n*n alkiota
        self.distance_matrix = [[0 for _ in range(n)] for _ in range(n)]
        self.propability_matrix = [[0 for _ in range(n)] for _ in range(n)]
        self.adjacency_matrix_weighted_cap = [[0 for _ in range(n)] for _ in range(n)]
        self.adjacency_matrix = [[0 for _ in range(n)] for _ in range(n)]
        self.capacity_matrix = [[0 for _ in range(n)] for _ in range(n)]
        self.min_cut_matrix = [[0 for _ in range(n)] for _ in range(n)]
        self.min_weight_degree = [[0 for _ in range(n)] for _ in range(n)]

        # create a Networkx graph from the nodes:
        self.G = nx.Graph()
        self.nodes = []
        i=0
        # next we create # n nodes with random coordinates on the square
        while i < self.n :
            noodi = (random.uniform(0, self.R), random.uniform(0, self.R))
            self.nodes.append(noodi)
            self.G.add_node(i, pos = noodi)
            i += 1


    def getDistance(self, node1: tuple, node2: tuple) -> float:
        # returns the distance between two nodes of the graph
        x = abs(node1[0] - node2[0])
        y = abs(node1[1] - node2[1])
        length = math.sqrt(x**2 + y**2)
        return length
    
    def distanceMatrix(self) -> list:
        #returns the distance matrix of the network
        for i in range(self.n):
            for j in range(self.n):
                self.distance_matrix[i][j] = self.getDistance(self.nodes[i], self.nodes[j])
        return self.distance_matrix

    def propabilityMatrix(self) -> list:
        # returns propability matrix from distance matrix according to p = e^(-D(x,y)/alfa*L)
        for i in range(self.n):
            for j in range(self.n):
                if i != j:
                    self.propability_matrix[i][j] = math.exp(- self.distance_matrix[i][j]/(self.alfa*self.R*math.sqrt(2)))
                else:
                    self.propability_matrix[i][j] = 0
        return self.propability_matrix
    
    def adjacencyMatrix(self) -> list:
        # determines links between nodes from the propabilities.
        # should the element be zero ore none when there is no link?
        for i in range(self.n):
            for j in range(i, self.n):
                if i != j:
                    k = random.uniform(0, 1.0)
                    if self.propability_matrix[i][j] >= k:
                        self.adjacency_matrix[i][j] = self.adjacency_matrix[j][i] = 1

                        distance = self.distance_matrix[i][j]
                        self.adjacency_matrix_weighted_cap[i][j] = self.adjacency_matrix_weighted_cap[j][i] = - math.log2(1 - 10**(-self.gamma*distance*0.1))

                        self.G.add_edge(i, j, weight=self.propability_matrix[i][j], capacity= self.adjacency_matrix_weighted_cap[i][j])
                    else:
                        self.adjacency_matrix[i][j] = self.adjacency_matrix[j][i] = None
                        self.adjacency_matrix_weighted_cap[i][j] = self.adjacency_matrix_weighted_cap[j][i] = None
                if i == j:
                    self.adjacency_matrix[i][j] = None
                    self.adjacency_matrix_weighted_cap[i][j] = None

        return self.adjacency_matrix
    
    def degree(self) -> list:
        # sum over rows gives degrees of nodes
        self.degrees = [sum(element for element in row if element is not None) for row in self.adjacency_matrix]
        return self.degrees
    
    def degree_weighted(self) -> list:
        # for convenience lets sum over rows so we get degrees of nodes weighted with capacity
        self.degree_weight = [sum(element for element in row if element is not None) for row in self.adjacency_matrix_weighted_cap]
        return self.degree_weight


    def capacityMatrix(self) -> list:
        #creates a capacity matrix from the distance_matrix according to maximum capacity of the link
        # C = log_2(1 - eta) = log_2(1-10^(-gamma*D/10))
        # oletuksena kaikki arvot matriisissa nollia, joten muutetaan arvoja vain kun noodien välillä linkki
        for i in range(self.n):
            for j in range(self.n):
                if self.adjacency_matrix[i][j] == 1:
                    distance = self.distance_matrix[i][j]
                    self.capacity_matrix[i][j] = - math.log2(1 - 10**(-self.gamma*distance*0.1))

        return self.capacity_matrix 
    
    def minimum_cuts(self) -> list:
        # palauttaa matriisin, jonka alkiot ovat minimileikkauksia noodien välillä
        for i in range(self.n):
            for j in range(i+1, self.n):
                if i != j:
                    self.min_cut_matrix[i][j] = self.min_cut_matrix[j][i] = nx.minimum_cut_value(self.G, i, j) 
        return self.min_cut_matrix
    
    def min_degree_weight(self) -> list:
        #palauttaa matriisin, jonka alkioina min(d_C(i),d_C(j)) missä d_C on noodin kapasiteetilla painotettu aste
        for i in range(self.n):
            for j in range(i+1, self.n):
                if i != j:
                    self.min_weight_degree[i][j] = self.min_weight_degree[j][i] = min(self.degree_weight[i],self.degree_weight[j])
        return self.min_weight_degree
    
