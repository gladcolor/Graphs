// CS786_A3_Huan.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

// list representation of graph 
#include <chrono>
#include <stdio.h> 
#include <stdlib.h> 
#include <limits.h> 
#include <vector>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include <stdio.h> 
#include <stdlib.h> 
#include <limits.h> 
#include <algorithm>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <math.h>       /* sqrt */
#include <random>

#include <list> 

using namespace std;
using namespace std::chrono;

struct edge {
	long iStart;
	long iEnd;
	double dW;
};

// Code modified from https://www.geeksforgeeks.org/breadth-first-search-or-bfs-for-a-graph/
// This class represents a directed graph using 
// adjacency list representation 
class Graph
{
	int V;    // No. of vertices 

	// Pointer to an array containing adjacency 
	// lists 
	list<int>* adj;

public:
	Graph(int V);  // Constructor 

	// function to add an edge to graph 
	void addEdge(int v, int w);

	// prints BFS traversal from a given source s 
	void BFS(int s);
};

Graph::Graph(int V)
{
	this->V = V;
	adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w); // Add w to v¡¯s list. 
}

void Graph::BFS(int s)
{
	// Mark all the vertices as not visited 
	bool* visited = new bool[V];
	for (int i = 0; i < V; i++)
		visited[i] = false;

	// Create a queue for BFS 
	list<int> queue;

	// Mark the current node as visited and enqueue it 
	visited[s] = true;
	queue.push_back(s);

	// 'i' will be used to get all adjacent 
	// vertices of a vertex 
	list<int>::iterator i;

	while (!queue.empty())
	{
		// Dequeue a vertex from queue and print it 
		s = queue.front();
		//cout << "Source: " << s << " ";
		queue.pop_front();

		// Get all adjacent vertices of the dequeued 
		// vertex s. If a adjacent has not been visited,  
		// then mark it visited and enqueue it 
		for (i = adj[s].begin(); i != adj[s].end(); ++i)
		{
			if (!visited[*i])
			{
				visited[*i] = true;
				queue.push_back(*i);
			}
		}
	}
}


int main(int argc, char** argv)
{
    std::cout << "Hello CS786!\n \n";
	 
	//ifstream read_file("sample_rmat_1m.gr");
	//ifstream read_file("sample_random100.gr");

	ifstream read_file(argv[1]);

	cout << "line:" << argv[1] << endl;

	string line;
	// skip the first 7 lines
	for (int i = 0; i < 7; i++)  
	{
		getline(read_file, line);
		//c
	}

	long n; 
	long m;
	long u;
	long v;
	double w;
	string t;  // tempary variable.
	read_file >> t;
	read_file >> t;

	read_file >> n;  
	read_file >> m;

	cout << "n:" << n << endl;
	cout << "m:" << m << endl;
		
	Graph g(n + 1);

	time_t start_time = time(NULL);

	clock_t ct;
	ct = clock();

	long lLine_cnt = 0;
	while (getline(read_file, line))
	{
		//cout << "line:" << line.c_str() << endl;
		read_file >> t;
		read_file >> u;
		read_file >> v;
		read_file >> w;

		g.addEdge(u, v);

		lLine_cnt++;

		if (lLine_cnt % 100000 == 0)
		{
			cout << "Current line: " << lLine_cnt << endl;
		}


		//cout << "u, v: " << u << ", " << v << endl;


	}

	time_t end_time = time(NULL);

	cout << "Seconds for data reading (s): " << (double)(clock() - ct) / CLOCKS_PER_SEC << endl;
	ct = clock();

	cout << "Constructed the graph."  << endl;
	cout << "Time used:  " << end_time - start_time << endl;

	start_time = time(NULL);


	for (int i = 0; i < 200; i++)
	{

		//int iDest = rand() % n;
		auto now = std::chrono::system_clock::now();
		auto now_ms = time_point_cast<std::chrono::milliseconds>(now);
		auto epoch = now_ms.time_since_epoch();
		auto value = duration_cast<std::chrono::milliseconds>(epoch);
		long duration = value.count();
		srand((unsigned int)duration);

		long value_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()).time_since_epoch()).count();

		std::random_device rd;
		std::mt19937 mt(rd());

		float a = (float)mt() / 2 / (float)INT_MAX;
		float b = (float)mt() / 2 / (float)INT_MAX;
		int iDest = a * n;
		int iSrc = b * n;

		printf("%d, iSrc: %d \n", i, iSrc);
		//int distance = dijkstra_2nodes(graph, iSrc, iDest);

		g.BFS(iSrc);

		now = std::chrono::system_clock::now();
		now_ms = time_point_cast<std::chrono::milliseconds>(now);
		epoch = now_ms.time_since_epoch();
		value = duration_cast<std::chrono::milliseconds>(epoch);
		duration = value.count();
		srand((unsigned int)(iDest * value_ms));

	}

 
	end_time = time(NULL);

	cout << "Time used:  " << end_time - start_time << endl;
	cout << "Seconds for BFS() 200 times (s): " << (double)(clock() - ct) / CLOCKS_PER_SEC << endl;

	cout << "Finished BFS." << endl;

	return 0;
	//for (int i = 0; i < 10; i++)
	//	in >> a[i];
	//for (int i = 0; i < 10; i++)
	//	cout << a[i] << endl;
 
}








// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
