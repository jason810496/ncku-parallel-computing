/*
Dijkstra with MPI

O(n^2) Dijkstra's algorithm with MPI
*/
#include<mpi.h>
#include<iostream>
#include<vector>
#include<utility>
#include<queue>

#define MAX_N 50005
#define F first
#define S second
#define INF 1e9

struct Edge{
    short t,w;
    Edge(short t,short w):t(t),w(w){}
    Edge(){}

    bool operator>(const Edge& e) const{
        return w > e.w;
    }
};


int main(int argc, char *argv[]){
    int n;
    int cluster_size;
    int worker_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    // MPI Data Type
    MPI_Datatype MPI_EDGE_TYPE;
    MPI_Type_contiguous(2,MPI_SHORT,&MPI_EDGE_TYPE);
    MPI_Type_commit(&MPI_EDGE_TYPE);

    vector< Edge > graph[MAX_N];
    vector< int > dis(MAX_N,INF);
    vector< bool > vis(MAX_N,false);
    priority_queue<Edge, vector<Edge>, greater<Edge>> pq;
    // input
    if(worker_id == 0){
        // file input
        char filename[100];
        scanf("%s",filename);
        FILE *fp = fopen(filename,"r");
        if (fp == NULL){
            printf("File not found\n");
            return 0;
        }
        // read n
        fscanf(fp,"%d",&n);
        for(int i=0;i<n;i++){
            short from,to,weight;
            fscanf(fp,"%hd %hd %hd",&from,&to,&weight);
            graph[from].push_back(Edge(to,weight));
        }
        fclose(fp);
    }
    // broadcast data
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(graph.data(),MAX_N,MPI_EDGE_TYPE,0,MPI_COMM_WORLD);
    // Dijkstra
    short start_vertex = 0;
    dis[start_vertex] = 0;

    while (!pq.empty()) {
        int u = pq.top().second; pq.pop();
        if (vis[u]) continue;
        vis[u] = true;
        // split all adj vertexes to workers
        // for (Edge e : adj[u]) {
        //     if (e.t < t || dist[e.v] <= dist[u] + e.w) continue;
        //     dist[e.v] = dist[u] + e.w;
        //     pq.push({dist[e.v], e.v});
        // }
    }


    MPI_Type_free(&MPI_EDGE_TYPE);
    MPI_Finalize();
    return 0;
}