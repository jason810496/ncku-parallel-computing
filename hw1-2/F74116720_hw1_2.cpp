#include<mpi.h>
#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<utility>

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


// debug print
void print_points(std::vector<Point> points){
    for(auto p: points){
        printf("(%d,%d) ",p.x,p.y);
    }
    printf("\n");
}

void print_ids(std::vector<Point> points){
    for(auto p: points){
        printf("%d ",p.id);
    }
    printf("\n");
}

#define COLLINEAR 0
#define CLOCKWISE 1
#define COUNTER_CLOCKWISE 2

// orientation of 3 points
u_int8_t orientation(Point p, Point q, Point r){
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if( val < 0) return COUNTER_CLOCKWISE;
    if ( val > 0) return CLOCKWISE;
    return COLLINEAR;
}

// compare 2 points with respect to pivot
bool compare(Point p1, Point p2, Point pivot){
    int o = orientation(pivot, p1, p2);
    if (o == COLLINEAR)
        return (std::sqrt((p1.x - pivot.x)*(p1.x - pivot.x) + (p1.y - pivot.y)*(p1.y - pivot.y)) <
                std::sqrt((p2.x - pivot.x)*(p2.x - pivot.x) + (p2.y - pivot.y)*(p2.y - pivot.y)));
    return (o == CLOCKWISE); // clockwise sort
}

// Helper function to compute the cross product of two vectors
int cross_product(const Point& p1, const Point& p2, const Point& p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

std::vector<Point> get_convex_hull(std::vector<Point> &points){
    // Andrew monotone chain
    int n = points.size();
    sort(points.begin(),points.end(),[](Point p1,Point p2){
        return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);
    });

    std::vector<Point> hull,lower_hull,upper_hull;
    
    for(int i=0;i<n;i++){
        while(lower_hull.size() >= 2 && cross_product(lower_hull[lower_hull.size()-2],lower_hull[lower_hull.size()-1],points[i]) <= 0){
            lower_hull.pop_back();
        }
        while (upper_hull.size() >= 2 && cross_product(upper_hull[upper_hull.size()-2],upper_hull[upper_hull.size()-1],points[i]) >= 0){
            upper_hull.pop_back();
        }
        lower_hull.push_back(points[i]);
        upper_hull.push_back(points[i]);
    }
    // merge the lower and upper hulls
    // hull = upper_hull from left to right + lower_hull from right to left
    for(int i=0;i<upper_hull.size();i++){
        hull.push_back(upper_hull[i]);
    }
    for(int i=lower_hull.size()-2;i>0;i--){ // don't include the first and last point
        hull.push_back(lower_hull[i]);
    }

    printf("get_convex_hull\n");
    print_ids(hull);
    return hull;
}


// Finds the index of the rightmost point in the hull
int find_rightmost_point(const std::vector<Point>& hull) {
    int idx = 0;
    for (int i = 1; i < hull.size(); i++) {
        if (hull[i].x > hull[idx].x || (hull[i].x == hull[idx].x && hull[i].y > hull[idx].y)) {
            idx = i;
        }
    }
    return idx;
}

// Finds the index of the leftmost point in the hull
int find_leftmost_point(const std::vector<Point>& hull) {
    int idx = 0;
    for (int i = 1; i < hull.size(); i++) {
        if (hull[i].x < hull[idx].x || (hull[i].x == hull[idx].x && hull[i].y < hull[idx].y)) {
            idx = i;
        }
    }
    return idx;
}

// Merge two convex hulls in a clockwise direction
std::vector<Point> merge_hulls(const std::vector<Point>& hull1, const std::vector<Point>& hull2) {
    std::vector<Point> merged_hull;

    // Get rightmost point of the first hull and leftmost point of the second hull
    int rightmost1 = find_rightmost_point(hull1);
    int leftmost2 = find_leftmost_point(hull2);

    // Add upper hull from hull1 (left to right)
    for (int i = 0; i <= rightmost1; ++i) {
        merged_hull.push_back(hull1[i]);
    }

    // Add upper hull from hull2 (left to right)
    for (int i = leftmost2; i < hull2.size(); ++i) {
        merged_hull.push_back(hull2[i]);
    }

    // Add lower hull from hull2 (right to left)
    for (int i = leftmost2 - 1; i >= 0; --i) {
        merged_hull.push_back(hull2[i]);
    }

    // Add lower hull from hull1 (right to left)
    for (int i = rightmost1 + 1; i < hull1.size(); ++i) {
        merged_hull.push_back(hull1[i]);
    }

    return merged_hull;
}

#define INT_MAX 2147483647
// add padding to the points
void add_padding(std::vector<Point>& points, int chunk_size){
    int n = points.size();
    int padding = chunk_size - n;
    while(padding--){
        points.push_back(Point(INT_MAX,INT_MAX,INT_MAX));
    }
}

// remove padding from the points
void remove_padding(std::vector<Point>& points){
    int idx = points.size()-1;
    while(points[idx].x == INT_MAX){
        points.pop_back();
        idx--;
    }
}


int main(int argc, char *argv[]){
    int n;
    int cluster_size;
    int worker_id;
    std::vector<Point> points(MAX_N);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &cluster_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);

    // MPI Data Type
    MPI_Datatype MPI_POINT_TYPE;
    MPI_Type_contiguous(3,MPI_INT,&MPI_POINT_TYPE);
    MPI_Type_commit(&MPI_POINT_TYPE);

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
            fscanf(fp,"%d %d",&points[i].x,&points[i].y);
            points[i].id = i+1;
        }
    }
    // broadcast data
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(points.data(),n,MPI_POINT_TYPE,0,MPI_COMM_WORLD);
    // divide the load
    int chunk_size = n / cluster_size;
    int start = worker_id * chunk_size;
    int end = (worker_id == cluster_size-1)? n: start + chunk_size;

    dbg(worker_id);
    dbg(start);
    dbg(end);
    dbg(points.size());
    // graham scan
    std::vector<Point> local_subset(points.begin()+start,points.begin()+end);
    std::vector<Point> local_hull = get_convex_hull(local_subset);
    // merge hulls
    std::vector<Point> global_hull;
    if(worker_id == 0){
        global_hull = local_hull;
        for(int i=1;i<cluster_size;i++){
            std::vector<Point> recv_hull(chunk_size);
            MPI_Recv(recv_hull.data(),chunk_size,MPI_POINT_TYPE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            print_ids(local_hull);
            remove_padding(recv_hull);

            dbg(worker_id);
            print_points(local_hull);
            print_points(global_hull);
            global_hull = merge_hulls(global_hull,recv_hull);
        }
    }
    else{
        if(local_hull.size() < chunk_size){
            add_padding(local_hull,chunk_size);
        }
        MPI_Send(local_hull.data(),chunk_size,MPI_POINT_TYPE,0,0,MPI_COMM_WORLD);
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