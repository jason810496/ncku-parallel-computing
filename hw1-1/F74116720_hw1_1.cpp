#include<iostream>
#include<cstdio>
#include<climits>
#include <mpi.h>

#define LL unsigned long long
#define DEBUG 1
#define dbg(v) std::cout << "Line(" << __LINE__ << ") -> " << #v << " = " << (v) << "\n";

void debug_bitmask(LL bitmask){
    if(DEBUG==0) return;
    for(int i=32;i>=0;i--){
        if(bitmask&(1<<i)){
            printf("1");
        }
        else{
            printf("0");
        }
    }
    printf("\n");
}

u_int16_t check(LL bitmask, u_int16_t n, u_int16_t m, u_int16_t val[], u_int16_t pos[][33]){
    LL mask = 0;
    // debug_bitmask(bitmask);
    for(u_int16_t i=0;i<n;i++){
        if(bitmask&(1<<i)){ // take this element
            u_int16_t j = 0;
            u_int16_t p;
            while(j <= m && pos[i][j]!= SHRT_MAX){
                p = pos[i][j];
                mask |= (1<<p);
                j++;
            }
        }
    }
    for(u_int16_t i=1;i<=n;i++){
        if(! (mask&(1<<i)) ){
            return 0;
        }
    }
    return 1;
}

int main (int argc, char *argv[]) {
    u_int16_t n,m;
    u_int16_t val[33]; // 1~32
    u_int16_t pos[33][33]; // 1~32
    int cluster_size;
    int worker_id;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);

    if(worker_id == 0){
        // file input
        char filename[100];
        scanf("%s",filename);
        FILE *fp = fopen(filename,"r");
        if (fp == NULL){
            printf("File not found\n");
            return 0;
        }
        fscanf(fp,"%hu %hu",&n,&m);
        // m input
        for(u_int16_t i=0;i<m;i++){
            u_int16_t cnt;
            fscanf(fp,"%hu %hu",&cnt,&val[i]);
            for(u_int16_t j=0;j<cnt;j++){
                fscanf(fp,"%hu",&pos[i][j]);
            }
            pos[i][cnt] = SHRT_MAX;
        }
        fclose(fp);
    }
    // worker 0 broadcast data
    MPI_Bcast(&n,1,MPI_UNSIGNED_SHORT,0,MPI_COMM_WORLD);
    MPI_Bcast(&m,1,MPI_UNSIGNED_SHORT,0,MPI_COMM_WORLD);
    MPI_Bcast(val,33,MPI_UNSIGNED_SHORT,0,MPI_COMM_WORLD);
    MPI_Bcast(pos,33*33,MPI_UNSIGNED_SHORT,0,MPI_COMM_WORLD);
    // all worker do the job
    LL total_size = 1ULL<<n;
    LL chunk_size = total_size / cluster_size;
    if (chunk_size == 0){
        chunk_size = 1;
    }
    LL current_bitmask = chunk_size * worker_id;
    LL end_bitmask;
    if(worker_id == cluster_size-1){
        end_bitmask = total_size;
    }
    else{
        end_bitmask = current_bitmask + chunk_size;
    }

    LL local_ans = 0;
    for(;current_bitmask<end_bitmask;current_bitmask++){
        local_ans += check(current_bitmask,n,m,val,pos);
    }
    LL global_ans = 0;
    // worker 0 gather data
    MPI_Reduce(&local_ans,&global_ans,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);

    if (worker_id == 0){
        printf("%llu\n",global_ans);
    }


    MPI_Finalize();
    return 0;
}