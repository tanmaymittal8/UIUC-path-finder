/*
  graph.h
  Project #7 - Openstreet Maps
  Author: Tanmay Mittal
  Date: 4/26/2022
  Class: CS 251, Spring 2022, University of Illinois at Chicago
  System: Replit on Windows
  Project discription: Basic graph class using adjacency matrix representation.  Currently limited to a graph with at most 100 vertices.
*/

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_map>
#include <map>
#include <algorithm>

using namespace std;

template<typename VertexT, typename WeightT>
class graph {
 private:
typedef unordered_map <VertexT, WeightT> adj_map;
set <VertexT> vertices;
unordered_map <VertexT, adj_map> adjList;

 public:

// default constructor, empty
  graph() {
  }

// returns the total # of vertices currently in graph
  int NumVertices() const {
    return vertices.size();
  }

// returns the total # of edges currently in graph
  int NumEdges() const {
    int count = 0;
    for (auto &e : adjList) {
      for (auto &i : e.second) {
        count += 1;
      }
    }
    return count;
  }

// adds a vertex v to graph if it dosen't exists
// returns false if vertex already exists
  bool addVertex(VertexT v) {
    if (adjList.count(v) > 0) {  // if v exists in graph
      return false;  // as already exists
    } else {  // add vertex
      vertices.emplace(v);
      unordered_map <VertexT, WeightT> temp;
        adjList.emplace(v, temp);
    }
    return true;
  }

// adds an edge from vertex "from" to vertex "to" and also
// updates the weight as well
// returns false if both vertices don't exist
  bool addEdge(VertexT from, VertexT to, WeightT weight) {
    if (adjList.count(from) == 0) {
      return false;  // if vertex from dosen't exists
    }
    if (adjList.count(to) == 0) {
      return false;  // if vertex to dosen't exists
    }
      for (auto &i : adjList.at(from)) {
    // find "to" vertex in the map
        if (i.first == to) {
          i.second = weight;  // set weight
          return true;
        }
      }
    adjList.at(from).emplace(to, weight);
    return true;
  }

// gets the weight of an edge. finds the edge from vertex
// "from" and vertex "to" and gets the weight of that edge
  bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
    if (adjList.count(from) == 0) {
      return false;  // if vertex from dosen't exists
    }
    if (adjList.count(to) == 0) {
      return false;  // if vertex to dosen't exists
    }
    // traverse the map
      for (auto &i : adjList.at(from)) {
        if (i.first == to) {
          weight = i.second;  // access the weight
          return true;
        }
      }
    return false;
  }

// gets all the neighbors of the parametrized vertex v
// finds all the vertices that have an edge from v
  set<VertexT> neighbors(VertexT v) const {
    set<VertexT> s;
    // if adjList map is empty for
    // vertex v so no edges from v
    if (adjList.count(v) == 0) {
      return s;
    }
    for (auto &e : adjList.at(v)) {
    // adds all the neighbours to a set
      s.emplace(e.first);
    }
    return s;
  }


  // ouptuts a vector of all the vertices available in the graph
  vector<VertexT> getVertices() const {
    vector <VertexT> node;
    // traverse the map
    for (auto &e : adjList) {
      node.push_back(e.first);
    }
    return node;  // returns a copy:
  }

  // dump
  // Dumps the internal state of the graph for debugging purposes.
  void dump(ostream& output) const {
    /*
    output << "***************************************************" << endl;
    output << "********************* GRAPH ***********************" << endl;

    output << "**Num vertices: " << this->NumVertices() << endl;
    output << "**Num edges: " << this->NumEdges() << endl;

    output << endl;
    output << "**Vertices:" << endl;
    for (int i = 0; i < this->NumVertices(); ++i) {
      output << " " << i << ". " << this->Vertices[i] << endl;
    }

    output << endl;
    output << "**Edges:" << endl;
    for (int row = 0; row < this->NumVertices(); ++row) {
      output << " row " << row << ": ";

      for (int col = 0; col < this->NumVertices(); ++col) {
        if (this->AdjMatrix[row][col].EdgeExists == false) {
          output << "F ";
        } else {
          output << "(T,"
            << this->AdjMatrix[row][col].Weight
            << ") ";
        }
      }
      output << endl;
    }
    output << "**************************************************" << endl;
    */
  }
};
