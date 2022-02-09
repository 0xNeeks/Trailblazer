/* Name: Nicholas Cuevas (SCPD Student)
 * Course: CS106B
 * Section Leader: Diwakar Ganesan (cganesan)
 * Assignment Description:
 * This assignment file is the cpp file of the Trailblazer Assignment.
 * In this file, you will find the function implementations for the Trailblazer class.
 * Additional Resources Used:
 * This Wikipedia page was used as additional context related to the shortest path problems. The page was used as a reference to understand
 * graph problems in greater depth, especially for finding the shortest path.
 * https://en.wikipedia.org/wiki/Shortest_path_problem#Strategic_shortest-paths
 * Additional Resources Used:
 * This Wikipedia page provided additional insights to the A* Algorithm.
 * https://en.wikipedia.org/wiki/A*_search_algorithm#Pseudocode
 * Additional Resources Used:
 * This research provided more background information and a general approach for
 * how to ignore edges for the alternate path algorithm.
 * https://www.informatik.hu-berlin.de/de/forschung/gebiete/wbi/research/publications/2015/sigspatial_kshortest.pdf
 * Additional Resources Used:
 * This Wikipedia page provided additional insights as to how to approach Kruskal's algorithm, specifically regarding clusters
 * https://en.wikipedia.org/wiki/Kruskal%27s_algorithm
 * Late Days Used: 1
*/

#include "trailblazer.h"
#include "queue.h"
#include "priorityqueue.h"

using namespace std;

/* Function Prototypes */
bool depthFirstSearchHelper(BasicGraph& graph, Vertex* start, Vertex* end, Vector<Vertex*>& path);
void constructPath(Map<Vertex*, Vertex*>& previousMap, Vector<Vertex*>& path, Vertex* end);
void DijkstraAStarAlternateHelper(BasicGraph& graph, Vertex* start, Vertex* end, Vector<Vertex*>& path, bool heuristic, Edge* edgeIgnored);
void mergeClusters(HashMap<Vertex*, HashSet<Vertex*> >& clusters, Vertex* start, Vertex* finish);

/*
 * The depthFirstSearch function passes a BasicGraph object by reference, and two Vertex pointers to start and end destinations.
 * The function resets the graph data. Then the function calls the recursive helper, setting the result (true if path found, false otherwise)
 * to the boolean pathFound. This boolean is used to determine whether to return an empty vector or the given path constructed.
 * Returns the path constructed through DFS search.
 */
Vector<Vertex*> depthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end) {
    graph.resetData();
    Vector<Vertex*> path;
    bool pathFound = depthFirstSearchHelper(graph, start, end, path);

    if (!pathFound) {
        path.clear();
    }

    return path;
}

/*
 * The depthFirstSearchHelper function passes a BasicGraph object by reference, two Vertex pointers to start and end destinations, and a Vector of vertices by reference (the path to be returned).
 * The function adds the start vertex to the path and sets its color to Green.
 * Checks if the start and end vertices are the same and returns true to path found when this is the case.
 * Otherwise a Set of Vertices is created from the current vertex's neighbors. The neighbors are traversed, checking if they have been visited.
 * If they have not been visited, the function recurses until it needs to unchoose and backtrack or has found a path.
 * Returns true if a path is found and false otherwise.
 */
bool depthFirstSearchHelper(BasicGraph& graph, Vertex* start, Vertex* end, Vector<Vertex*>& path) {
    path.add(start); //choose
    start->setColor(GREEN);

    if (start == end) {
        return true;
    }

    Set<Vertex*> neighborsOfVertex = graph.getNeighbors(start);
    for (Vertex* neighbor : neighborsOfVertex) {
        if (neighbor->getColor() == UNCOLORED) { //if the neighbor has not been visited
            if (depthFirstSearchHelper(graph, neighbor, end, path)) { //explore
                return true;
            }
        }
    }

    path.removeBack(); //unchoose
    start->setColor(GRAY);
    return false;
}

/*
 * The breadthFirstSearch function passes a BasicGraph object by reference, and two Vertex pointers to start and end destinations.
 * The function resets the graph data. Then the function adds the start vertex to the queue of visitedVertices and sets its color to Yellow.
 * Checks if the start and end vertices are the same and adds the end location to the path, then breaks the loop.
 * Otherwise checks the neighbors and if they haven't been visited yet, enqueues the vertex, setting the color to yellow and updating the previous
 * vertex map.
 * Reconstructs the path based on the map.
 * Returns the path constructed through BFS search.
 */
Vector<Vertex*> breadthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end) {
    graph.resetData();

    Vector<Vertex*> path;
    Map<Vertex*, Vertex*> previousMap; //for tracking the path
    Queue<Vertex*> visitedVerteces;

    visitedVerteces.enqueue(start); //initial vertex in the queue
    start->setColor(YELLOW);

    while (!visitedVerteces.isEmpty()) {
        Vertex* current = visitedVerteces.dequeue();
        current->setColor(GREEN);

        if (current == end) {
            path.add(end);
            break;
        }

        Set<Vertex*> neighborsOfVertex = graph.getNeighbors(current);

        //checks the neighboring vertices, enqueuing unvisited ones
        for (Vertex* neighbor : neighborsOfVertex) {
            if (neighbor->getColor() != YELLOW && neighbor->getColor() != GREEN) { //if the neighbor has not been visited
                visitedVerteces.enqueue(neighbor);
                neighbor->setColor(YELLOW);
                previousMap[neighbor] = current; //tracks the previous of a vertex
            }
        }
    }

    constructPath(previousMap, path, end); //reconstruct the path based on the map
    path.insert(0, start);

    return path;
}

/*
 * The constructPath function passes a Map object of Vertex* as the key and value by reference, a Vector of vertices path by reference, and two Vertex pointers to start and end destinations.
 * The function begins at the end location of the path. Checks if that vertex as a key is in the map and that it was visited in the search. The value of
 * that key is inserted at the front of the path vector.
 * The function then recurses, using the current key as the end vertex.
 */
void constructPath(Map<Vertex*, Vertex*>& previousMap, Vector<Vertex*>& path, Vertex* end) {
    Vertex* current = previousMap.get(end);

    if (previousMap.containsKey(current) && current->getColor() == GREEN) {
        path.insert(0,current);
        constructPath(previousMap, path, current);
    }
}

/*
 * The dijkstrasAlgorithm function passes a BasicGraph object by reference, and two Vertex pointers to start and end destinations.
 * The function resets the graph data. Then the function sets the heuristic boolean to false
 * indicating to the helper function that it does not use the heuristic function.
 * Also passes nullptr to the helper indicating no edges are to be ignored.
 * The helper function for Dijkstra, A*, and Alternate path Algorithms is then called to create the path.
 * Returns the path constructed through Dijkstras algorithm.
 */
Vector<Vertex*> dijkstrasAlgorithm(BasicGraph& graph, Vertex* start, Vertex* end) {
    graph.resetData();
    Vector<Vertex*> path;
    bool heuristic = false;

    DijkstraAStarAlternateHelper(graph, start, end, path, heuristic, nullptr);

    return path;
}

/*
 * The aStar function passes a BasicGraph object by reference, and two Vertex pointers to start and end destinations.
 * The function resets the graph data. Then the function sets the heuristic boolean to true
 * indicating to the helper function that it uses the heuristic function.
 * Also passes nullptr to the helper indicating no edges are to be ignored.
 * The helper function for Dijkstra, A*, and Alternate path Algorithms is then called to create the path.
 * Returns the path constructed through A* algorithm.
 */
Vector<Vertex*> aStar(BasicGraph& graph, Vertex* start, Vertex* end) {
    graph.resetData();
    Vector<Vertex*> path;
    bool heuristic = true;

    DijkstraAStarAlternateHelper(graph, start, end, path, heuristic, nullptr);

    return path;
}

/*
 * The DijkstraAStarAlternateHelper function passes a BasicGraph object by reference, two Vertex pointers to start and end destinations, a boolean indicating whether its a heuristic,
 * and an Edge pointer to the edge to be ignored.
 * The function heuristic will be true for AStar and Alternate path algorithms call this helper, but false when used for Dijkstra's algorithm. The Edge pointer for the edge to be ignored
 * will be nullptr for Dijkstra and AStar algorithms. Will pass the edge to ignore for alternate path.
 * The function creates a HashMap to track the vertex priority of the vertices in the graph. This is updated as the algorithm searches.
 * The function adds the start vertex with its start cost to the prioirty queue of pqueue and sets its color to Yellow.
 * Dequeue from the prioirty queue, set color to Green to indicate visited the vertex and then
 * checks if the start and end vertices are the same and adds the end location to the path, then breaks the loop.
 * Then while the priority queue is not empty, checks the current vertex neighbors. If they haven't been visited before,
 * enqueues the neighbor vertex if the total cost of the current edge is less than the neighbor's total cost (initially postiive infinity).
 * Adds the neighbor vertex to the prioirty queue with its prioirty if not already in the queue, otherwise changes its prioirty.
 * Also checks if its not an edge to ignore, which the conditional will always be true since nullptr is passed. False is only possible
 * for the altenrate path algorithm when an edge to ignore is passed.
 * Reconstructs the path based on the map.
 * Returns the path constructed through the search.
 * Assumes valid parameters are passed for the appropriate algorithm use.
 */
void DijkstraAStarAlternateHelper(BasicGraph& graph, Vertex* start, Vertex* end, Vector<Vertex*>& path, bool heuristic, Edge* edgeIgnored) {
    Map<Vertex*, Vertex*> previousMap;
    PriorityQueue<Vertex*> pqueue;
    HashSet<Vertex*> vertexInpQueue;
    HashMap<Vertex*, double> vertexPriority;

    //adds vertices with a positive infinity priority to the hashmap
    for (Vertex* vertex : graph) {
        vertexPriority.add(vertex, POSITIVE_INFINITY);
    }

    double cost = 0.0; //cost for Dijkstra

    //cost for A* and alternate path algorithms
    if (heuristic) {
        cost = cost + heuristicFunction(start, end);
    }

    //initial vertex
    pqueue.enqueue(start, cost);
    vertexInpQueue.add(start);
    start->setColor(YELLOW);

    while (!pqueue.isEmpty()) {
        double currentPriority = pqueue.peekPriority();
        Vertex* current = pqueue.dequeue();
        vertexInpQueue.remove(current);
        current->setColor(GREEN);

        if (current == end) {
            path.add(end);
            break;
        }

        //checks the neighboring vertices, enqueuing unvisited ones
        Set<Vertex*> neighborsOfVertex = graph.getNeighbors(current);
        for (Vertex* neighbor : neighborsOfVertex) {

            if (neighbor->getColor() != YELLOW && neighbor->getColor() != GREEN) { //if the neighbor has not been visited
                Edge* currentEdge = graph.getEdge(current, neighbor);
                double edgeCost = currentEdge->cost;
                double totalCost = edgeCost + currentPriority;

                if (totalCost < vertexPriority[neighbor]) {
                    //Heuristic for new cost for A* and Altenrate path algorithms
                    if (heuristic) {
                        totalCost = totalCost + heuristicFunction(neighbor, end);
                    }

                    if (vertexInpQueue.contains(neighbor) && (currentEdge != edgeIgnored)) { //edge ignored conditional always true unless called by alternate path, then possibly false
                       pqueue.changePriority(neighbor, totalCost); //changes priority of neighbor vertex in the priority queue with appropriate totalCost
                       vertexPriority[neighbor] = totalCost; //updates priority in the hashmap tracking vertex priorities
                       previousMap[neighbor] = current; //tracks the previous of a vertex

                    } else if (!vertexInpQueue.contains(neighbor) && (currentEdge != edgeIgnored)) {
                        pqueue.enqueue(neighbor, totalCost); //adds neighbor to the priority queue with appropriate totalCost
                        vertexInpQueue.add(neighbor); //tracks if vertex already in the prioirty queue
                        vertexPriority[neighbor] = totalCost; //updates priority in the hashmap tracking vertex priorities
                        previousMap[neighbor] = current; //tracks the previous of a vertex
                        neighbor->setColor(YELLOW);
                    }
                }
            }
        }
    }

    constructPath(previousMap, path, end); //reconstruct the path based on the map
    path.insert(0, start);

}

/*
 * The alternatePath function passes a BasicGraph object by reference,  two Vertex pointers to start and end destinations, and a double representing the difference.
 * The function resets the graph data. Then the function calls the algorithm helper function with the parameters for it to be an A* shortest path search.
 * The helper function call updates the shortestPathVector.
 * Reset the graph data to find the alternate path.
 * The vertices of the shortest path are stored in a hashset to simplify the difference calculation later in the algorithm
 * For each vector in the shortest path, the edge to be ignored is found. Then the algorithm helper function is called, indicating
 * that the heuristic is true (so an A* call) but also indicates to ignore the given edge.
 * The found alternate path has its vertices stored in a HashSet to calculate the difference.
 * The calculated difference must be larger than the given difference in order to add it to the prioirty queue.
 * The priority queue stores the alternate path's cost as that path's priority.
 * If the shortest path isn't empty, the highest priority path is the next best alternate path that is statistically different.
 * That path is set to the path of vertices to be returned by the function.
 * Assumes valid passed parameters..
 * Returns the path constructed through the Alternate Path algorithm.
 */
Vector<Vertex*> alternatePath(BasicGraph& graph, Vertex* start, Vertex* end, double difference) {
    graph.resetData();

    Vector<Vertex*> path;
    PriorityQueue<Vector<Vertex*> > shortestAltPath;
    Vector<Vertex*> shortestPathVector;

    //find the shortest path
    DijkstraAStarAlternateHelper(graph, start, end, shortestPathVector, true, nullptr); //uses the AStar algorithm based on its parameters

    graph.resetData(); //reset the graph without the AStar data

    //store the shortest path vertices in a set to simplify the difference calculation later
    HashSet<Vertex*> shortestPathSet;
    for (Vertex* vertex : shortestPathVector) {
        shortestPathSet.add(vertex);
    }

    //For each edge in the best path it will produce one candidate alternate route
    for (int i=0; i < shortestPathVector.size()-1; i++) {
        Edge* edgeRemoved = graph.getEdge(shortestPathVector[i], shortestPathVector[i+1]);
        Vector<Vertex*> altPathVector;
        altPathVector.clear(); //clear vector for each edge iteration

        DijkstraAStarAlternateHelper(graph, start, end, altPathVector, true, edgeRemoved); //helper function called for each edge in the best path, ignoring one at a time

        HashSet<Vertex*> altPathSet;
        for (Vertex* vertex : altPathVector) {
            altPathSet.add(vertex);
        }

        //calculate the difference between the best and alternate paths
        double nodesInAltNotInBestPath = (double) (altPathSet - shortestPathSet).size();
        double altDifference = nodesInAltNotInBestPath / (double) altPathSet.size();

        if (altDifference > difference) {
            double altShortestPathCost = 0.0;

            //find the alternate path cost
            for (int i=0; i < altPathVector.size()-1; i++) {
                Edge* edge = graph.getEdge(altPathVector[i], altPathVector[i+1]);
                altShortestPathCost += edge->cost;
            }

            shortestAltPath.enqueue(altPathVector, altShortestPathCost); //enqueue the alternate path with its cost as its priority
        }
    }

    //set the highest priority, alternate path to the vector to be returned if the prioirty queue is not empty
    if (!shortestAltPath.isEmpty()) {
      path = shortestAltPath.dequeue();
    }

    return path;
}

/*
 * The kruskal function passes a BasicGraph object by reference.
 * The function resets the graph data. Then the function creates a set of all the graph's edges.
 * These edges are then added to the priority queue based on their cost.
 * A HashMap is created mapping a vertex, to a hashset of all the vertices connected. This is intiialized with the vertices connected to itself.
 * While the prioirty queue has 2 or more clusters, the current edge is dequeued and its start and finish vertices are set.
 * The function then checks the hashmap if the start vertex does not contain the finish vertex in the value-hashSet.
 * If this is the case, the clusters are merged between the start and finished vertices and the edge is added to the Set of edges
 * to be returned.
 * The mergeClusters helper function merges them.
 * Returns the set of edges to create the mst.
 * Assumption and Note: Works well for small and medium sized files, but takes too long for large files (though it does create an mst)
 */
Set<Edge*> kruskal(BasicGraph& graph) {
    graph.resetData();

    Set<Edge*> mst;
    PriorityQueue<Edge*> edgesQueue;
    Set<Edge*> tempEdgesSet = graph.getEdgeSet();

    //initialize priority queue
    for (Edge* edge : tempEdgesSet) {
        double edgeCost = edge->cost;
        edgesQueue.enqueue(edge, edgeCost);
    }

    HashMap<Vertex*, HashSet<Vertex*> > clusters;
    int clusterCount = 0;

    //initialize clusters
    for (Vertex* vertex : graph) {
        clusters[vertex].add(vertex);
        clusterCount++;
    }

    while (clusterCount >= 2) {
        Edge* currentEdge = edgesQueue.dequeue(); //work with the current edge
        Vertex* start = currentEdge->start;
        Vertex* finish = currentEdge->finish;

        //start vertex does not contain the finish vertex in its cluster
        if (!clusters[start].contains(finish)) {
            mergeClusters(clusters, start, finish);
            mst.add(currentEdge);
            clusterCount--;
        }
    }

    return mst;
}

/*
 * The mergeClusters function passes a hashmap of a vertex pointer key and hashset of vertex pointers key by reference. Also passes a Vertex pointer for the start and finish.
 * The function combines the hashsets of the start and finish clusters.
 * Then for each vertex in the combinedCluster, checks if the clusters hashmap contains that vertex as a key. If so, sets that vertex to the comvbinedCluster
 * Assumption and Note: Works well for small and medium sized files, but takes too long for large files (though it does create an mst)
 */
void mergeClusters(HashMap<Vertex*, HashSet<Vertex*> >& clusters, Vertex* start, Vertex* finish) {
    HashSet<Vertex*> combinedCluster = clusters[start] + clusters[finish];

    for (Vertex* vertex : combinedCluster) {
        if (clusters.containsKey(vertex)) {
            clusters[vertex] = combinedCluster;
        }
    }
}
