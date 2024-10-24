/*
Dijkstra with MPI

O(n^2) Dijkstra's algorithm with MPI
*/
#include<mpi.h>
#include<iostream>
#include<vector>
#include<utility>
#include<queue>
#include<string.h>

#define MAX_N 50001
#define MAX_N_2 50001L*50001L
#define F first
#define S second
#define INF 1e9

#define min(a,b) ((a)<(b)?(a):(b))

void dijksra(int start_node,int end_node,short graph[MAX_N][MAX_N],int local_dis[MAX_N]){
    std::priority_queue<std::pair<int,int>,std::vector<std::pair<int,int>>,std::greater<std::pair<int,int>>> pq;
    pq.push({0,start_node});
    while(!pq.empty()){
        auto [d,u] = pq.top();
        pq.pop();
        if(d > local_dis[u]) continue;
        for(int v=0;v<MAX_N;v++){
            if(graph[u][v] == 0) continue;
            if(local_dis[v] > local_dis[u] + graph[u][v]){
                local_dis[v] = local_dis[u] + graph[u][v];
                pq.push({local_dis[v],v});
            }
        }
    }
}


int main(int argc, char *argv[]){
    int n;
    int cluster_size;
    int worker_id;
    short graph[MAX_N][MAX_N];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);

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
            graph[from][to] = weight;
        }
        fclose(fp);
        // init data
        memset(graph,0,sizeof(graph));
    }
    // broadcast data
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    // MPI_Bcast(graph,MAX_N_2,MPI_SHORT,0,MPI_COMM_WORLD);
    // divide the load
    int chunk_size = n / cluster_size;
    int start = worker_id * chunk_size;
    int end = (worker_id == cluster_size-1)? n: start + chunk_size;

    // init local_dis for all workers
    int local_dis[MAX_N];
    for(int i=0;i<n;i++){
        local_dis[i] = INF;
    }
    local_dis[0] = 0;
    dijksra(start,end,graph,local_dis);
    // global result
    int global_dis[MAX_N];

    // Dijkstra
    if( worker_id == 0 ){
        // init global_dis
        for(int i=0;i<n;i++){
            global_dis[i] = INF;
        }
        global_dis[0] = 0;
        // self result
        for(int i=0;i<n;i++){
            global_dis[i] = min(global_dis[i],local_dis[i]);
        }
        // merge the results 
        for(int i=1;i<cluster_size;i++){
            short recv_dis[MAX_N];
            MPI_Recv(recv_dis,MAX_N,MPI_SHORT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            for(int j=0;j<n;j++){
                global_dis[j] = min(global_dis[j],recv_dis[j]);
            }
        }
    }
    else{
        // send the results to the worker 0
        MPI_Send(local_dis,MAX_N,MPI_SHORT,0,0,MPI_COMM_WORLD);
    }

    // print the result
    if(worker_id == 0){
        for(int i=0;i<n;i++){
            printf("%d ",global_dis[i]);
        }
    }

    return 0;
}