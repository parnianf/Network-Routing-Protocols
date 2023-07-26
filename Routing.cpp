#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <set>
#include <chrono>
using namespace std;


const int MAXN = 100;
const float INFMAX = 100000;
const float INFMIN = -1e9;

const string TOPOLOGY_CMD = "topology";
const string SHOW_CMD = "show";
const string LSRP_CMD = "lsrp";
const string DVRP_CMD = "dvrp";
const string MODIFY_CMD = "modify";
const string REMOVE_CMD = "remove";

struct Edge {
	int u, v;
	bool operator ==(const Edge& e) {
		return (u == e.u) && (v == e.v);
	}

};

typedef vector<Edge> Edges;

class LinkState {

private:
	float dist[MAXN];
	int parent[MAXN];

public:
	int getMinNodeInDijkstra(bool visited[], int numberOfNodes) {
		int node = 1, tempMin = INFMAX;
		for (int i = 1; i <= numberOfNodes; i++) {
			if ((!visited[i]) && tempMin >= dist[i]) {
				node = i;
				tempMin = dist[i];
			}
		}
		return node;
	}

	void dijkstra(int src, int numberOfNodes, float cost[MAXN][MAXN]) {
		for (int i = 0; i < MAXN; i++) {
			dist[i] = INFMAX;
			parent[i] = -1;
		}

		dist[src] = 0;

		bool visited[MAXN] = { 0 };
		int iterNum = 1;
		for (int i = 1; i <= numberOfNodes; i++) {
			if (i == src)
				continue;
			int nearest = getMinNodeInDijkstra(visited, numberOfNodes);
			visited[nearest] = 1;
			for (int adj = 1; adj <= numberOfNodes; adj++) {
				float newDistance = dist[nearest] + cost[nearest][adj];
				if (!visited[adj] && dist[nearest] != INFMAX && cost[nearest][adj] != -1 && newDistance < dist[adj]) {
					dist[adj] = newDistance;
					parent[adj] = nearest;
				}
			}
			printIterations(iterNum, numberOfNodes);
			iterNum++;
		}
		print(numberOfNodes, src);
	}

	void printIterations(int iterNum, int numberOfNodes) {
		cout << "Iter: " << iterNum << endl;
		cout << "Dest |";
		for (int n = 1; n <= numberOfNodes; n++) {
			cout << std::setw(3) << n << " |";
		}
		cout << endl;
		cout << "Cost |";
		for (int n = 1; n <= numberOfNodes; n++) {
			if (dist[n] == INFMAX)
				cout << std::setw(3) << "-1" << " |";
			else
				cout << std::setw(3) << dist[n] << " |";
		}
		cout << endl;
		for (int i = 0; i < numberOfNodes; i++) cout << "--------";
		cout << endl;
	}


	void lsrp(int src, int numberOfNodes, float cost[MAXN][MAXN], Edges edges) {
		set <int> vertices;
		for (auto edge : edges)
			vertices.insert(edge.u);


		if (src != -1) {
			auto start = std::chrono::steady_clock::now();
			dijkstra(src, numberOfNodes, cost);
			auto end = std::chrono::steady_clock::now();
			std::cout << "Convergence time of LSRP: " << std::chrono::duration <double, std::milli>(end - start).count() << " ms" << std::endl;

		}
		else {
			for (auto v : vertices) {
				dijkstra(v, numberOfNodes, cost);
			}
		}
	}

	void printPath(int v, string& path) {
		if (parent[v] == -1) {
			path += to_string(v) + "->";
			return;
		}
		printPath(parent[v], path);
		path += to_string(v) + "->";
	}

	void print(int numberOfNodes, int src) {
		cout << "Path: [s] -> [d] | \tMin-Cost \t|\t Shortest Path";
		for (int i = 1; i < numberOfNodes + 1; i++) {
			if (i == src)
				continue;

			cout << endl;
			string path = "";
			printPath(i, path);
			path.pop_back();
			path.pop_back();
			cout << "\t" << src << " -> " << i << " \t\t " << dist[i] << " \t\t\t " << path;
		}
		cout << endl;
		for (int i = 0; i < numberOfNodes; i++) cout << "-------------------";
		cout << endl;
	}

};


class DistanceVector {

private:
	float dist[MAXN];
	int parent[MAXN];
	int nextHop;
	bool nextHopFound;

public:
	DistanceVector() {
		nextHopFound = false;
	}

	void print(int src, int numberOfNodes)
	{
		cout << "Dest	\tNext Hop	Dist	\tShortest Path" << endl;
		for (int i = 0; i < numberOfNodes; i++) cout << "-------------------";
		cout << endl;
		for (int i = 1; i <= numberOfNodes; i++) {
			string path = "";
			nextHop = src;
			printPath(i, path);
			path.pop_back();
			path.pop_back();
			cout << i << "\t\t" << nextHop << "\t\t" << dist[i] << "\t\t[" << path << "]" << endl;
		}
		cout << endl;
	}

	void printPath(int v, string& path) {
		if (parent[v] == -1) {
			nextHopFound = true;
			path += to_string(v) + "->";
			return;
		}
		printPath(parent[v], path);
		if (nextHopFound) {
			nextHop = v;
			nextHopFound = false;
		}
		path += to_string(v) + "->";
	}


	void bellmanFord(int src, int numberOfNodes, vector<Edge> edges, float cost[MAXN][MAXN])
	{
		for (int i = 0; i < MAXN; i++) {
			dist[i] = INFMAX;
			parent[i] = -1;
		}
		dist[src] = 0;

		for (int i = 1; i <= numberOfNodes - 1; i++)
		{
			for (int j = 0; j < edges.size(); j++)
			{
				int u = edges[j].u;
				int v = edges[j].v;
				int weight = cost[u][v];
				if (dist[u] != INFMAX && (dist[u] + weight < dist[v])) {
					dist[v] = dist[u] + weight;
					parent[v] = u;
				}
			}
		}
		print(src, numberOfNodes);
	}

	void dvrp(int src, int numberOfNodes, vector<Edge> edges, float cost[MAXN][MAXN]) {
		set <int> vertices;
		for (auto edge : edges)
			vertices.insert(edge.u);

		if (src != -1) {
			auto start = std::chrono::steady_clock::now();
			bellmanFord(src, numberOfNodes, edges, cost);
			auto end = std::chrono::steady_clock::now();
			std::cout << "Convergence time of DVRP: " << std::chrono::duration <double, std::milli>(end - start).count() << " ms" << std::endl;
		}
		else {
			for (auto v : vertices) {
				bellmanFord(v, numberOfNodes, edges, cost);
			}
		}
	}
};


class Routing {
public:

	void initialize() {
		for (int i = 0; i < MAXN; i++) {
			for (int j = 0; j < MAXN; j++) {
				if (i == j)
					cost[i][j] = 0;
				else
					cost[i][j] = -1;
			}
		}
		edges.clear();
		numberOfNodes = 0;
	}

	vector<string> parse(const string& input, char delim) {
		vector<string> parsed;
		stringstream ssInp(input);
		string str;
		while (getline(ssInp, str, delim))
			parsed.push_back(str);
		return parsed;
	}

	void setTopology(const vector<string>& args) {
		numberOfNodes = INFMIN;
		for (int i = 0; i < args.size(); i++) {
			vector<string> edgeTpg = parse(args[i], '-');
			edges.push_back({ stoi(edgeTpg[0]), stoi(edgeTpg[1]) });
			edges.push_back({ stoi(edgeTpg[1]), stoi(edgeTpg[0]) });
			int u = stoi(edgeTpg[0]);
			int v = stoi(edgeTpg[1]);
			int maximum = max(u, v);
			numberOfNodes = max(maximum, numberOfNodes);
			cost[u][v] = stof(edgeTpg[2]);
			cost[v][u] = stof(edgeTpg[2]);
		}
	}

	void showTopology() {
		cout << " u|v  ";
		for (int i = 1; i <= numberOfNodes; i++) {
			cout << std::setw(3) << i;
		}
		cout << endl;

		for (int i = 0; i<numberOfNodes; i++) {
			cout << std::setw(3) << "------";
		}
		cout << endl;
		for (int i = 1; i <= numberOfNodes; i++) {
			cout << std::setw(3) << i << " | ";
			for (int j = 1; j <= numberOfNodes; j++) {
				cout << std::setw(3) << cost[i][j];
			}
			cout << endl;
		}
	}

	void modifyTopology(string topologies_str) {
		vector<string> topology_str = parse(topologies_str, '-');
		int source = stoi(topology_str[0]);
		int dest = stoi(topology_str[1]);
		float weight = stof(topology_str[2]);
		if (dest == source) {
			cout << "ERROR: Source and destination cannot be the same!" << endl;
			return;
		}
		numberOfNodes = max(numberOfNodes, max(source, dest));
		if (cost[source][dest] == -1) {
			cost[source][dest] = cost[dest][source] = weight;
			cost[source][source] = cost[dest][dest] = 0;
			edges.push_back({ source, dest });
			edges.push_back({ dest, source });
		}
		else {
			cost[source][dest] = cost[dest][source] = weight;
		}
	}

	void removeEdge(int src, int dest) {
		for (int i = 0; i < edges.size(); i++) {
			if (edges[i].u == src && edges[i].v == dest) {
				edges.erase(std::remove(edges.begin(), edges.end(), edges[i]), edges.end());
			}
		}
	}

	void removeTopology(string removeCommand) {
		vector<string> parsedCommand = parse(removeCommand, '-');
		int source = stoi(parsedCommand[0]);
		int dest = stoi(parsedCommand[1]);

		if (dest == source) {
			cout << "ERROR: Source and destination cannot be the same!" << endl;
			return;
		}

		if (cost[source][dest] == -1) {
			cout << "ERROR: this topology does NOT exists!" << endl;
			return;
		}

		cost[source][dest] = cost[dest][source] = -1;
		removeEdge(source, dest);
		removeEdge(dest, source);
	}

	void handleCommand() {
		string input;
		while (getline(cin, input)) {
			vector<string> parsedCommand = parse(input, ' ');
			string command = parsedCommand[0];
			vector<string> arguments(parsedCommand.begin() + 1, parsedCommand.end());

			if (command == TOPOLOGY_CMD) {
				initialize();
				setTopology(arguments);
			}

			else if (command == SHOW_CMD) {
				showTopology();
			}

			else if (command == LSRP_CMD) {
				LinkState* ls = new LinkState();
				if (arguments.size() == 1) {
					ls->lsrp(stoi(arguments[0]), numberOfNodes, cost, edges);
				}
				else {
					ls->lsrp(-1, numberOfNodes, cost, edges);
				}
			}
			else if (command == DVRP_CMD) {
				DistanceVector* dv = new DistanceVector();
				if (arguments.size() == 1)
					dv->dvrp(stoi(arguments[0]), numberOfNodes, edges, cost);
				else
					dv->dvrp(-1, numberOfNodes, edges, cost);
			}
			else if (command == MODIFY_CMD) {
				modifyTopology(arguments[0]);
			}
			else if (command == REMOVE_CMD) {
				removeTopology(arguments[0]);
			}
		}
	}

private:
	int numberOfNodes;
	Edges edges;
	float cost[MAXN][MAXN];
};


int main() {
	Routing routing;
	routing.handleCommand();
	return 0;
}
