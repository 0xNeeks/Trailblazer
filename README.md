# Trailblazer 

Project is part of Stanford's CS 106B: Programming Abstractions Course and the code authored can be found in src/trailblazer.cpp

This assignment focuses on graphs, specifically on searching for paths in a graph.
- depth­first search (DFS) 
- breadth­first search (BFS) 
- Dijkstra's algorithm
- A* search
- Alternate Path

Problem Description:
This program displays various 2­dimensional worlds that represent either maps, mazes, or terrain and allows 
the user to generate paths in a world from one point to another. When you start up the program, you will see a graphical window containing a 2D maze, 
where white squares are open and black ones represent walls. The program is also able to display terrain, where bright colors indicate higher elevations 
and darker colors represent lower elevations. Mountain ranges appear in bright white, while deep canyons are closer to black.

If you click on any two points in the world, the program will find a path from the starting position to the ending position. As it does so, it will color 
the vertexes green, yellow, and gray based on the colors assigned to them by the algorithm. Once the path is found, the program will highlight it and display 
information about the path weight in the console. The user can select one of five path­searching algorithms in the top menu:
- depth­first search (DFS) 
- breadth­first search (BFS) 
- Dijkstra's algorithm
- A* search
- Alternate Path
