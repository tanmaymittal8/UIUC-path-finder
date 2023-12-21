/*
  application.cpp
  Project - Openstreet Maps
  Author: Tanmay Mittal
  Date: 4/26/2022
  Project discription: The file traverses the graphs and finds the closest meeting plce of two individuals and gives the total distance required to travel and the path with node ID's.
  Creative component: Finds the distance and outputs a path for a person to go from one place to another. It also has a functunality of adding a new stop on the trip and finding the total distance and path to reach both the destenations. The user should choose option c when prompted and can input any graph provided. Then the user is required to enter the starting and ending destination. After the path and distance is outputted the user can enter another destination to add on to the original destination and then the program will output the total distance and the path from the last destination to the new destination. To exit press #.
*/

// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//


#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <stack>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include "graph.h"
#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"

using namespace std;
using namespace tinyxml2;
const double INF = numeric_limits<double>::max();  // global variable

class prioritize {
 public:
  bool operator()(
    const pair<long long, double>& p1,
    const pair<long long, double>& p2) const{
    return p1.second > p2.second;
  }
};

// set up a priority queue
priority_queue<
  pair<long long, double>,          // (key,value) pair
  vector<pair<long long, double>>,  // store pairs in a vector
  prioritize> unvisitedQueue;   // function object


// adds all the buildings in the graph into a set
void add_buildings_set(
    set<string> &buildings_left, vector<BuildingInfo>& Buildings) {
  buildings_left.clear();  // resets the vector back to 0
  for (auto &y : Buildings) {  // adds all the buildings
      buildings_left.emplace(y.Fullname);
    }
}

// search for buildings from user input
BuildingInfo searchBuilding(vector<BuildingInfo> Buildings,
                            string personBuilding, bool &found) {
  for (auto p : Buildings) {
         if (personBuilding == p.Abbrev) {  // input mathches abbreb
           found = true;
           return p;  // return building
         } else if (p.Fullname.find(personBuilding) != string::npos) {
           found = true;
           return p;  // return building
    }  // checks, input matches fullname or part of it
  }
    found = false;  // if building not found
    return {};
}

// find the nearist building betw persn 1 and persn 2
BuildingInfo nearestBuilding(
    vector<BuildingInfo> Buildings, Coordinates midpoint,
    set<string> buildings_left) {
  double min = INF;
  BuildingInfo buildingCenter;
  for (auto &e : Buildings) {
    // checks if this building has already been ckecked and no path found
    if (buildings_left.count(e.Fullname) != 0) {
        double distance = distBetween2Points(midpoint.Lat,
                      midpoint.Lon, e.Coords.Lat, e.Coords.Lon);
        if (distance < min) {  // min algorithm
          min = distance;
          buildingCenter = e;  // save the min dist building
        }
      }
    }
  cout << " " << buildingCenter.Fullname << endl;
  cout << " (" << buildingCenter.Coords.Lat << ", ";
  cout << buildingCenter.Coords.Lon << ")" << endl;
  return buildingCenter;
}



// find nearest node to a building
long long nearestNode(BuildingInfo b, vector<FootwayInfo>& Footways,
                map<long long, Coordinates>& Nodes) {
  double min = INF;
  long long nearestNode;
  for (auto &e : Footways) {
     for (int i = 0; i < e.Nodes.size(); i++) {  // each node in footway
      double lat1 = Nodes.at(e.Nodes.at(i)).Lat;
      double lng1 = Nodes.at(e.Nodes.at(i)).Lon;
      double distance = distBetween2Points(lat1, lng1, b.Coords.Lat, b.Coords.Lon);
      if (distance < min) {  // min algorithm
         min = distance;
      nearestNode = Nodes.at(e.Nodes.at(i)).ID;  // save nearest node
     }
    }
  }
  cout << " " << nearestNode << endl;
  cout << " (" << Nodes.at(nearestNode).Lat << ", ";
  cout << Nodes.at(nearestNode).Lon << ")" << endl;
  return nearestNode;
}

// get the path from start to destination
vector<long long> getPath(map<long long, long long> predecessors,
                    long long nearestNode2) {
  long long currV;
  stack<long long> path;
  vector <long long> pathp;
  currV = nearestNode2;  // traverse backwards, dest
  while (currV != 0) {
    path.push(currV);
    // traverse the predecessor map until it reached 0
    currV = predecessors.at(currV);
  }
  // inverse direction of the input
  while (path.size() != 0) {
    currV = path.top();
    path.pop();
    pathp.push_back(currV);
  }
  return pathp;
}

// Dijkstra's algorithm
void DijkstraShortestPath(graph<long long, double> G, long long nearestNode1, map<long long, long long> &predecessors, map<long long, double> &distances) {
    set <long long > visited1;
    vector<long long> allNodes = G.getVertices();
    for (auto &vertex : allNodes) {  // for each through every vertex
      distances[vertex] = INF;  // initial distance to infinity
      predecessors[vertex] = 0;  // initialize
      unvisitedQueue.push(make_pair(vertex, INF));  // push node to queue
    }
    distances[nearestNode1] = 0;  // initialize
    unvisitedQueue.push(make_pair(nearestNode1, 0));
    while (unvisitedQueue.size() != 0) {  // not empty
      pair<long long, double> currNode = unvisitedQueue.top();
      unvisitedQueue.pop();
      if (currNode.second == INF) {
      break;  // can't visit currntV
    } else if (visited1.count(currNode.first) != 0) {
        continue;  // currentV has been visited
    }  // visit current v
      visited1.insert(currNode.first);
       set<long long> neighbors = G.neighbors(currNode.first);
    for (auto &neighbor : neighbors) {
      double edgeWeight = 0;
      G.getWeight(currNode.first, neighbor, edgeWeight);
      double pathDist = currNode.second + edgeWeight;
      if (pathDist < distances[neighbor]) {
        distances[neighbor] = pathDist;
        predecessors[neighbor] = currNode.first;
        unvisitedQueue.push(make_pair(neighbor, pathDist));
      }  // set distance from start vertex to end for each vertex
    }  // predecessor array track path from each vertex
      }
}

// prints and calculates the output path in terminal for 
// both people
void print_path(
  map<long long, double> distances1, map<long long, double> distances2,
  vector <long long> printPath1, vector <long long> printPath2, map<long long,
  long long> predecessors1, map<long long, long long> predecessors2,
  long long nearestNode1Center) {
  cout << "Person 1's distance to dest: " << distances1.at(nearestNode1Center);
  cout << " miles" <<  endl;
  printPath1 = getPath(predecessors1, nearestNode1Center);
  cout << "Path: ";
  for (int i = 0; i < printPath1.size()-1; i++) {
    cout << printPath1.at(i) << "->";

  }  // the last path Id won't have an arrow
  cout << printPath1.at(printPath1.size()-1) << endl;
  cout << "Person 2's distance to dest: " << distances2.at(nearestNode1Center);
  cout << " miles" <<  endl;
  printPath2 = getPath(predecessors2, nearestNode1Center);
  cout << "Path: ";
  for (int c = 0; c < printPath2.size()-1; c++) {
    cout << printPath2.at(c) << "->";
  }
  cout << printPath2.at(printPath2.size()-1) << endl << endl;
}

// prints the intitial position and buildings of both individuals
void print_initial_locat(BuildingInfo b1, BuildingInfo b2) {
  cout << "Person 1's point:" << endl;
        cout << " " << b1.Fullname << " " <<  endl;
        cout << " (" << b1.Coords.Lat << ", " << b1.Coords.Lon << ")"<< endl;
        cout << "Person 2's point:" << endl;
        cout << " " << b2.Fullname << " " <<  endl;
        cout << " (" << b2.Coords.Lat << ", " << b2.Coords.Lon << ")"<< endl;
}

// traverses the graph and finds the destination building and all the nodes associated with it
// finds the distance and path to the most optimal shortes distance possible 
// or finds the next best destination untill none of the buildins can be reached so it 
// gives an message that the destination isn't possible
void traverse_graph(bool &found, BuildingInfo b1, BuildingInfo b2,
      BuildingInfo buildingCenter, map<long long, Coordinates>& Nodes,
      vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings,
      string &person2Building, Coordinates &midpoint, long long nearestNode1,
      long long nearestNode2, long long nearestNode1Center,
      map<long long, long long> predecessors1, map<long long, long long>
      predecessors2, map<long long, double> distances1,
      map<long long,double> distances2, graph<long long, double> G,
      vector <long long> printPath1, vector <long long> printPath2,
      set <string> buildings_left, string filename) {
  if (found) {  // checks if b1 found in graph, or error statnm
    b2 = searchBuilding(Buildings, person2Building, found);
      if (found) {  // checks if b2 found in graph, or error statnm

          print_initial_locat(b1, b2);
          // locate center building
          midpoint = centerBetween2Points(b1.Coords.Lat, b1.Coords.Lon,
                                        b2.Coords.Lat, b2.Coords.Lon);
          cout << "Destination Building: " << endl;
          buildingCenter = nearestBuilding(Buildings, midpoint, buildings_left);
          cout << endl;
          // FINDS NEAREST NODES FROM BUILDINGS 1,2 & CENTER
          
          cout << "Nearest P1 node:" << endl;
          nearestNode1 = nearestNode(b1,  Footways, Nodes);
          cout << "Nearest P2 node:" << endl;
          nearestNode2 = nearestNode(b2,  Footways, Nodes);
          cout << "Nearest destination node:" << endl;
          nearestNode1Center = nearestNode(buildingCenter,  Footways, Nodes);
      // RUN DIJKSTRA'S ALGORITHM
          DijkstraShortestPath(G, nearestNode1, predecessors1, distances1);
          DijkstraShortestPath(G, nearestNode2, predecessors2, distances2);
          if (distances1[nearestNode2] >= INF) {  // if no path possible
              cout << "Sorry, destination unreachable. " << endl;
            } else if (distances1[nearestNode1Center] >= INF || distances2[nearestNode1Center] >= INF) {  // find a different next closest path until no path is possibe
                while (buildings_left.size() == 0 || (distances1[nearestNode1Center] >= INF
                      && distances2[nearestNode1Center] >= INF)) {
            cout << "At least one person was unable to ";
            cout << "reach the destination building. ";
            cout << "Finding next closest building..." << endl;
            // keeps track of all the buildings tested
            buildings_left.erase(buildingCenter.Fullname);
            cout << "New destination building:" << endl;
            // recalculate center building
            buildingCenter = nearestBuilding(Buildings, midpoint, buildings_left);
            cout << "Nearest destination node:" << endl;
            nearestNode1Center = nearestNode(buildingCenter,  Footways, Nodes);
          }
            print_path(distances1, distances2, printPath1, printPath2,
                  predecessors1, predecessors2, nearestNode1Center);
              } else {  // print path
              print_path(distances1, distances2, printPath1, printPath2,
                  predecessors1, predecessors2, nearestNode1Center);
              }
        } else { cout << "Person 2's building not found" << endl; }
      } else { cout << "Person 1's building not found" << endl; }
}

// calculates and outputs the total distance and path to the destination
void application(map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings, graph<long long, double> G, string filename) {
  vector <long long> printPath1, printPath2;
  set <string> buildings_left;
  add_buildings_set(buildings_left, Buildings);
  string person1Building, person2Building;
  cout << endl;
  for (const auto& element : buildings_left) {
      cout << element << ", ";
  }
  cout << endl << endl;
  cout << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  while (person1Building != "#") {
    cout << "Enter person 2's building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    BuildingInfo b1, b2, buildingCenter;
    bool found = false;
    long long nearestNode1, nearestNode2, nearestNode1Center;
    map<long long, long long> predecessors1, predecessors2;
    map<long long, double> distances1, distances2;
    b1 = searchBuilding(Buildings, person1Building, found);
    Coordinates midpoint;
    traverse_graph(found, b1, b2, buildingCenter, Nodes, Footways,
      Buildings, person2Building, midpoint, nearestNode1, nearestNode2,
      nearestNode1Center, predecessors1, predecessors2, distances1,
      distances2, G, printPath1, printPath2, buildings_left, filename);
    add_buildings_set(buildings_left, Buildings);  // reset buildings_left
    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);
      }
}

// find the path from one building to another with all the
// buildings close to the path mentioned to eaiser navigation
void creative(
    map<long long, Coordinates>& Nodes, vector<FootwayInfo>& Footways,
    vector<BuildingInfo>& Buildings, graph<long long, double> G) {
  string person1Building, person2Building;
  set <string> buildings_left;
  add_buildings_set(buildings_left, Buildings);
    for (const auto& element : buildings_left) {
      cout << element << ", ";
  }
  cout << endl << endl;
  cout << endl << "Enter person 1's building (partial name or abbreviation), or #> ";
  getline(cin, person1Building);
  long long nearestNode1, nearestNode2;
  map<long long, long long> predecessors1;
  map<long long, double> distances1;
  double total_distance = 0;
  while (person1Building != "#") {
    cout << endl << endl << "Enter next building (partial name or abbreviation)> ";
    getline(cin, person2Building);
    vector <long long> printPath;
    BuildingInfo b1, b2;
    bool found = false;
    b1 = searchBuilding(Buildings, person1Building, found);
    b2 = searchBuilding(Buildings, person2Building, found);
    nearestNode1 = nearestNode(b1,  Footways, Nodes);
    nearestNode2 = nearestNode(b2,  Footways, Nodes);
    DijkstraShortestPath(G, nearestNode1, predecessors1, distances1);
    printPath = getPath(predecessors1, nearestNode2);
    total_distance =  total_distance + distances1.at(nearestNode2);
    cout << "Person 1's distance to dest: " << total_distance << endl;
    for (int i = 0; i < printPath.size()-1; i++) {
      cout << printPath.at(i) << "-> ";
    }
      cout << printPath.at(printPath.size()-1);
      person1Building = person2Building;
    }
}

int main() {
  // maps a Node ID to it's coordinates (lat, lon)
  map<long long, Coordinates>  Nodes;
  // info about each footway, in no particular order
  vector<FootwayInfo>          Footways;
  // info about each building, in no particular order
  vector<BuildingInfo>         Buildings;
  XMLDocument                  xmldoc;

  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  string def_filename = "map.osm";
  string filename;
  cout << "Choose from file name uic.osm or uiuc.osm" << endl;
  cout << "Enter map filename> ";
  getline(cin, filename);

  if (filename == "") {
    filename = def_filename;
  }

  //
  // Load XML-based map file
  //
  if (!LoadOpenStreetMap(filename, xmldoc)) {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return 0;
  }

  //
  // Read the nodes, which are the various known positions on the map:
  //
  int nodeCount = ReadMapNodes(xmldoc, Nodes);

  //
  // Read the footways, which are the walking paths:
  //
  int footwayCount = ReadFootways(xmldoc, Footways);

  //
  // Read the university buildings:
  //
  int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);

  //
  // Stats
  //
  assert(nodeCount == (int)Nodes.size());
  assert(footwayCount == (int)Footways.size());
  assert(buildingCount == (int)Buildings.size());

  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;


// adds vertices
// loops through the Nodes and adds all the verties in it
  graph <long long, double> G;
  for (auto &e : Nodes) {
    G.addVertex(e.first);
  }

// adds edges
  for (auto &t : Footways) {  // loops through footways
    // finds the nodes in each footway and finds the distance
    // btw two nodes
    for (long long i = 0; i < t.Nodes.size()-1; i++) {
      double lat1 = Nodes.at(t.Nodes.at(i)).Lat;
      double lng1 = Nodes.at(t.Nodes.at(i)).Lon;
      double lat2 = Nodes.at(t.Nodes.at(i+1)).Lat;
      double lng2 = Nodes.at(t.Nodes.at(i+1)).Lon;
      double distance = distBetween2Points(lat1, lng1, lat2, lng2);
      // adds the edges in both directions 
      G.addEdge(t.Nodes.at(i), t.Nodes.at(i+1), distance);
      G.addEdge(t.Nodes.at(i+1), t.Nodes.at(i), distance);
    }
  }
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;

  // Menu
  string userInput;
  cout << "Enter \"a\" for the standard application or "
        << "\"c\" for the creative component application> ";
  getline(cin, userInput);
  if (userInput == "a") {
    application(Nodes, Footways, Buildings, G, filename);
  } else if (userInput == "c") {
  // creative component
    creative(Nodes, Footways, Buildings, G);
  }
  cout << "** Done **" << endl;
  return 0;
}
