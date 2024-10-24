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

#define MAX_N 1500L
#define MAX_N_2 MAX_N*MAX_N
#define F first
#define S second
#define INF 1e9

#define min(a,b) ((a)<(b)?(a):(b))

short graph[MAX_N][MAX_N];
int dis[MAX_N];

struct Message{
    int start_node;
    int end_node;
    int from;

    Message(int start_node,int end_node,int from):start_node(start_node),end_node(end_node),from(from){}
    Message(){}
};

Message create_message(int n,int worker_size,int worker_id,int from){
    int chunk_size = n / worker_size;
    int start = worker_id * chunk_size;
    int end = (worker_id == worker_size-1 ? n: start + chunk_size);
    printf("create msg: start = %d, end = %d, from = %d\n",start,end,from);
    return Message(start,end,from);
}

Message create_exit_message(){
    return Message(-1,-1,-1);
}

int main(int argc, char *argv[]){
    int n;
    int cluster_size;
    int worker_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);

    MPI_Datatype MPI_MESSAGE_TYPE;
    MPI_Type_contiguous(2,MPI_INT,&MPI_MESSAGE_TYPE);
    MPI_Type_commit(&MPI_MESSAGE_TYPE);

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
        while(!feof(fp)){
            short from,to,weight;
            fscanf(fp,"%hd %hd %hd\n",&from,&to,&weight);
            printf("from = %d, to = %d, weight = %d\n",from,to,weight);
            graph[from][to] = weight;
        }
        fclose(fp);
        // init dis
        for(int i=0;i<n;i++){
            dis[i] = INF;
        }
        dis[0] = 0;
    }
    // broadcast data
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(graph,MAX_N_2,MPI_SHORT,0,MPI_COMM_WORLD);
    MPI_Bcast(dis,MAX_N,MPI_INT,0,MPI_COMM_WORLD);
    // divide the load
    int worker_size = cluster_size - 1;
    int chunk_size = n / worker_size;

    // init dis
    if( worker_id == 0){
        // dijkstra
        std::priority_queue<std::pair<int,int>,std::vector<std::pair<int,int> >,std::greater<std::pair<int,int> > > pq;
        pq.push({0,0});
        while(!pq.empty()){
            std::pair<int,int> pii = pq.top();
            pq.pop();
            int u = pii.S;
            int d = pii.F;
            printf("pop: u = %d, d = %d, dis[u] = %d\n",u,d,dis[u]);
            if(d > dis[u]) continue;
            // distribute the work to other workers
            /*original*/
            /*
            for(int v=0;v<MAX_N;v++){
                if(graph[u][v] == 0) continue;
                if(local_dis[v] > local_dis[u] + graph[u][v]){
                    local_dis[v] = local_dis[u] + graph[u][v];
                    pq.push({local_dis[v],v});
                }
            }
            */

            
            std::vector<int> local_dis(MAX_N);
            for(int i=0;i<n;i++){
                local_dis[i] = dis[i];
            }
            if( n > worker_size){
                std::vector<Message> msgs(worker_size);
                for(int ith_worker=0;ith_worker<worker_size;ith_worker++){
                        msgs[ith_worker].from = u;
                        msgs[ith_worker].start_node = ith_worker * chunk_size;
                        msgs[ith_worker].end_node = (ith_worker == worker_size-1 ? n: msgs[ith_worker].start_node + chunk_size);
                        MPI_Send(&msgs[ith_worker],1,MPI_MESSAGE_TYPE,ith_worker+1,0,MPI_COMM_WORLD);
                        MPI_Send(local_dis.data(),MAX_N,MPI_INT,ith_worker+1,0,MPI_COMM_WORLD);
                }

                std::vector<int> recv_dis(MAX_N);
                for(int ith_worker=0;ith_worker<worker_size;ith_worker++){
                    printf("coordinator waiting %d\n",ith_worker+1);
                    MPI_Recv(recv_dis.data(),MAX_N,MPI_INT,ith_worker+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    printf("coordinator recv from %d\n",ith_worker+1);
                    for(int v=msgs[ith_worker].start_node;v<msgs[ith_worker].end_node;v++){
                        // printf("v = %d, recv_dis = %d, dis = %d\n",v,recv_dis[v],dis[v]);
                        if(dis[v] > recv_dis[v]){
                            dis[v] = recv_dis[v];
                            printf("coordinator: push %d %d\n",dis[v],v);
                            pq.push({dis[v],v});
                        }
                    }
                }
                printf("coordinator: pq size = %ld\n",pq.size());
            
                // update dis
                // MPI_Bcast(dis,MAX_N,MPI_INT,0,MPI_COMM_WORLD);

                printf("coordinator: broadcast done\n");
            }
            else{ // no need to distribute the work
                printf("coordinator: no need to distribute the work\n");
                for(int v=0;v<n;v++){
                    printf("graph[%d][%d] = %d\n",u,v,graph[u][v]);
                    if(graph[u][v] == 0) continue;
                    printf("u = %d, v = %d, dis[v] = %d, dis[u] = %d, wt = %d\n",u,v,dis[v],dis[u],graph[u][v]);
                    if(dis[v] > dis[u] + (int)graph[u][v]){
                        dis[v] = dis[u] + (int)graph[u][v];
                        pq.push({dis[v],v});
                        printf("coordinator: push v:%d, dis[%d] = %d\n",v,v,dis[v]);
                    }
                }
            }

            // MPI_Bcast(dis,MAX_N,MPI_INT,0,MPI_COMM_WORLD);
        }
        printf("coordinator: dijkstra done\n");

        // send exit message after the dijkstra is done
        for(int ith_worker=0;ith_worker<worker_size;ith_worker++){
            Message msg = create_exit_message();
            MPI_Send(&msg,1,MPI_MESSAGE_TYPE,ith_worker+1,0,MPI_COMM_WORLD);
        }

        printf("coordinator: exit message sent\n");
    }
    else{ // workers
        while(true){
            std::vector<int> local_min(MAX_N);
            // for(int i=0;i<MAX_N;i++){
            //     local_min[i] = dis[i];
            // }
            Message msg;

            printf("worker %d: waiting\n",worker_id);
            MPI_Recv(&msg,1,MPI_MESSAGE_TYPE,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(msg.start_node == -1 || msg.end_node == -1 || msg.from == -1){
                break;
            }

            MPI_Recv(local_min.data(),MAX_N,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            printf("worker %d: start = %d, end = %d, from = %d\n",worker_id,msg.start_node,msg.end_node,msg.from);
            

            int from = msg.from;
            for(int v=msg.start_node;v<msg.end_node;v++){
                if(graph[from][v] == 0) continue;
                printf("worker %d: from = %d, v = %d, local_min = %d, graph = %d, local_min + graph = %d\n",worker_id,from,v,local_min[from],graph[from][v],local_min[from] + graph[from][v]);
                if(local_min[v] > local_min[from] + graph[from][v]){
                    printf("worker %d find %d %d %d\n",worker_id,local_min[from],graph[from][v],local_min[from] + graph[from][v]);
                    local_min[v] = local_min[from] + graph[from][v];
                }
            }
            MPI_Send(local_min.data(),MAX_N,MPI_INT,0,0,MPI_COMM_WORLD);
            printf("worker %d: send to coordinator\n",worker_id);

        }
        printf("worker %d: exit\n",worker_id);
    }

    // print the result
    if(worker_id == 0){
        printf("dis\n");
        for(int i=0;i<n;i++){
            printf("%d ",dis[i]);
        }
    }

    MPI_Type_free(&MPI_MESSAGE_TYPE);
    MPI_Finalize();

    return 0;
}