#include<mpi.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

#define MAX_N 15000
#define LL unsigned long long
#define DEBUG 1
#define dbg(v) std::cout << "Line(" << __LINE__ << ") -> " << #v << " = " << v << std::endl;

/*
    convex hull using graham scan with MPI
*/

struct Point{
    int id,x,y;
    Point(int id,int x,int y):id(id),x(x),y(y){}
    Point(){}
};

#define COLLINEAR 0
#define CLOCKWISE 1
#define COUNTER_CLOCKWISE 2

// orientation of 3 points
u_int8_t orientation(Point p, Point q, Point r){
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if(val == 0) return COLLINEAR;
    return (val > 0)? CLOCKWISE: COUNTER_CLOCKWISE;
}

// compare 2 points with respect to pivot
bool compare(Point p1, Point p2, Point pivot){
    u_int8_t o = orientation(pivot,p1,p2);
    if(o == COLLINEAR){
        return (std::sqrt(std::pow(pivot.x-p1.x,2) + std::pow(pivot.y-p1.y,2)) >= std::sqrt(std::pow(pivot.x-p2.x,2) + std::pow(pivot.y-p2.y,2)));
    }
    return o == CLOCKWISE;
}

// graham scan
std::vector<Point> graham_scan(std::vector<Point> points){
    int n = points.size();
    if(n < 3) return {};

    // find the bottommost point
    Point pivot = *std::min_element(points.begin(),points.end(),[](Point p1, Point p2){
        return (p1.y < p2.y) || (p1.y == p2.y && p1.x < p2.x);
    });

    // sort the points based on polar angle
    std::sort(points.begin(),points.end(),[pivot](Point p1, Point p2){
        return compare(p1,p2,pivot);
    });

    std::vector<Point> hull;
    hull.push_back(points[0]);
    hull.push_back(points[1]);
    hull.push_back(points[2]);

    for(int i=3;i<n;i++){
        while(orientation(hull[hull.size()-2],hull[hull.size()-1],points[i]) != COUNTER_CLOCKWISE){
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }

    return hull;
}

// merge 2 convex hulls
std::vector<Point> merge_hulls(std::vector<Point>& hull1, std::vector<Point>& hull2) {
    // Combine the points from both hulls and compute the convex hull of the union
    std::vector<Point> all_points = hull1;
    all_points.insert(all_points.end(), hull2.begin(), hull2.end());
    return graham_scan(all_points);
}



int main(int argc, char *argv[]){
    int n;
    int cluster_size;
    int worker_id;
    std::vector<Point> points;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);

    // MPI Data Type
    MPI_Datatype MPI_POINT;
    MPI_Type_contiguous(3,MPI_INT,&MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

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
        points.resize(n);
        for(int i=0;i<n;i++){
            fscanf(fp,"%d %d",&points[i].x,&points[i].y);
            points[i].id = i+1;
        }
    }
    // broadcast data
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(points.data(),n,MPI_POINT,0,MPI_COMM_WORLD);
    // divide the load
    int chunk_size = n / cluster_size;
    int start = worker_id * chunk_size;
    int end = (worker_id == cluster_size-1)? n: start + chunk_size;

    dbg(worker_id);
    dbg(start);
    dbg(end);
    dbg(points.size());
    // graham scan
    std::vector<Point> local_hull = graham_scan(std::vector<Point>(points.begin()+start,points.begin()+end));
    // merge hulls
    std::vector<Point> global_hull;
    if(worker_id == 0){
        global_hull = local_hull;
        for(int i=1;i<cluster_size;i++){
            std::vector<Point> recv_hull(chunk_size);
            MPI_Recv(recv_hull.data(),chunk_size*sizeof(Point),MPI_BYTE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            global_hull = merge_hulls(global_hull,recv_hull);
        }
    }
    else{
        MPI_Send(local_hull.data(),local_hull.size()*sizeof(Point),MPI_BYTE,0,0,MPI_COMM_WORLD);
    }
    // print the result
    if(worker_id == 0){
        for(auto p: global_hull){
            printf("%d ",p.id);
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}