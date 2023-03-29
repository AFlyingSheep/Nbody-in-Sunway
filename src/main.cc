#include <cstdio>
#include <math.h>
#include <cstdlib>
//#include <mpi.h>
#include <assert.h>
#include <cmath>
#include <vector>
#include <algorithm>
// #include <crts.h>

#include "main.h"
#include "utils.h"
#include "zone.h"

const real_type G = 6.67e-5;

const real_type timestep = 0.1;
const int cycles_times = 2;
const int ITA = 0.001;
#define FILE_NAME "particles2.txt"

real_type delta_t = 1 / timestep;
// Zone* z;

typedef struct {
    real_type px[1024];
    real_type py[1024];
    real_type pz[1024];
    real_type vx[1024];
    real_type vy[1024];
    real_type vz[1024];

    real_type m[1024];
} sendbuf;

const int nei_dx[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
const int nei_dy[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

void calculate_a(Body& ba, Body& bb) {

    real_type dx = ba.px_ - bb.px_, dy = ba.py_ - bb.py_;
    // printf("dxdydz%lf%lf%lf", dx, dy, dz);
    real_type d_2 = dx * dx + dy * dy;
    real_type d = sqrt(d_2);

    // printf("%lf %lf %lf", d, d_2, dx);

    // printf("dxdydz%lf", d);
    // F = GMm / r3
    // Fx = GMm / r2 * (rx / r)
    // ax = GM / r2 * (rx / r)

    ba.ax_ -= G * bb.m_ / (pow(d, 3) + ITA) * dx;
    ba.ay_ -= G * bb.m_ / (pow(d, 3) + ITA) * dy;
}

void calculate_a_zone(Zone* z) {
    for (size_t i = 0; i < z->size_; i++) {
        z->bodys_[i].ax_ = 0;
        z->bodys_[i].ay_ = 0;
        // z->bodys_[i].az_ = 0;
        for (size_t j = 0; j < z->size_; j++) {
            if (i == j) continue;
            calculate_a(z->bodys_[i], z->bodys_[j]);
        }
    }
}

void calculate_a_zone(Zone* z, size_t start, size_t end) {
    // assert(start <= end);
    for (size_t i = start; i < end; i++) {
        z->bodys_[i].ax_ = 0;
        z->bodys_[i].ay_ = 0;
        // z->bodys_[i].az_ = 0;
        for (size_t j = 0; j < z->size_; j++) {
            if (i == j) continue;
            calculate_a(z->bodys_[i], z->bodys_[j]);
        }
    }
}

inline void calculate_v(Body& b, real_type max_v) {
    b.vx_ += b.ax_ * delta_t;
    b.vy_ += b.ay_ * delta_t;
    // b.vz_ += b.az_ * delta_t;
    if (std::abs(b.vx_) > max_v) b.vx_ = (b.vx_ > 0)? max_v : -1.0 * max_v;
    if (std::abs(b.vy_) > max_v) b.vy_ = (b.vy_ > 0)? max_v : -1.0 * max_v;
}

inline void calculate_v_zone(Zone* z, real_type max_v) {
    for (size_t i = 0; i < z->size_; ++i) {
        calculate_v(z->bodys_[i], max_v);
    }
}

inline void calculate_v_zone(Zone* z, size_t start, size_t end, real_type max_v) {
    for (size_t i = start; i < end; ++i) {
        calculate_v(z->bodys_[i], max_v);
    }
}

inline void calculate_p(Body& b) {
    b.px_ += b.vx_ * delta_t;
    b.py_ += b.vy_ * delta_t;
    // b.pz_ += b.vz_ * delta_t;
}

inline void calculate_p_zone(Zone* z, std::vector<Body>& out_bodys) {
    for (size_t i = 0; i < z->size_; i++) {
        Body& b = z->bodys_[i];
        b.px_ += b.vx_ * delta_t;
        b.py_ += b.vy_ * delta_t;
        
        if (b.px_ < z->zone_low.x || b.py_ < z->zone_low.y || 
            b.px_ > z->zone_high.x || b.py_ > z->zone_high.y) {
            out_bodys.push_back(b);
            z->pop(i);
        }
    }
}

void calculate_test(const Zone* z, size_t start, size_t end, real_type dt) {
    for (size_t i = start; i < end; i++) {
        z->bodys_[i].px_ += dt;
        z->bodys_[i].py_ += dt;
        // z->bodys_[i].pz_ += dt;
    }
}

void print_body(Zone* z,bool is_all, size_t start = 0, size_t end = 0) {
    if (is_all) {
        start = 0;
        end = z->size_;
    }
    for (int i = start; i < end; i++) {
        z->bodys_[i].print_result();
    }
}

bool check(Zone* z1, Zone* z2) {
    // printf("z1: %d, z2: %d\n", z1->size_, z2->size_);
    assert(z1->size_ == z2->size_);

    // print_body(z2, true);
    for (size_t i = 0; i < z1->size_; i++) {
        Body& b1 = z1->bodys_[i], &b2 = z2->bodys_[i];
        if (b1 == b2) continue;
        return false;
    }
    return true;
}

template <typename T>
struct Pos;

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);
    // CRTS_init();

    double timeSt;

    // getchar();
    // FILE *file = fopen("output.txt", "w");
    srand(time(NULL));

    // input data
    FILE *file = fopen(FILE_NAME, "r");
    if (file == nullptr) {
        printf("Open file error!\n");
        return 1;
    }

    real_type x, y, vx, vy;
    real_type l0;
    int num_particles;

    fscanf(file, "%d %lf", &num_particles, &l0);
    Body* bodies = (Body*)malloc(sizeof(Body) * num_particles);

    real_type max_v = l0 / delta_t;

    Pos<real_type> pos_max, pos_min;
    pos_max.x = -1000;
    pos_max.y = -1000;
    pos_min.x = 1000;
    pos_min.y = 1000;

    for (int i = 0; i < num_particles; i++) {
        fscanf(file, "%lf %lf %lf %lf", &x, &y, &vx, &vy);
        // printf("粒子 %d 的坐标为 (%lf, %lf)\n", i+1, x, y);
        bodies[i] = Body(x, y, vx, vy);
        // bodies[i].print_result();

        pos_max.x = std::max(pos_max.x, x);
        pos_max.y = std::max(pos_max.y, y);

        pos_min.x = std::min(pos_min.x, x);
        pos_min.y = std::min(pos_min.y, y);
    }

    pos_min.x -= l0 / 2;
    pos_min.y -= l0 / 2;
    pos_max.x += l0 / 2;
    pos_max.y += l0 / 2;

    real_type dx = pos_max.x - pos_min.x;
    real_type dy = pos_max.y - pos_min.y;

    // printf("dx: %lf, dy: %lf\n", dx, dy);

    // num of zones
    int zone_num_x = std::ceil(dx / l0);
    int zone_num_y = std::ceil(dy / l0);

    // **判断越界
    int zone_num = zone_num_x * zone_num_y;

    Zone* zones = (Zone*)malloc(sizeof(Zone) * zone_num_x * zone_num_y);

    for (int i = 0; i < zone_num; i++) {
        zones[i] = Zone();
        Pos<int> p = h2d(i, zone_num_x, zone_num_y);
        zones[i].zone_low.x = p.x * l0 + pos_min.x;
        zones[i].zone_low.y = p.y * l0 + pos_min.y;

        zones[i].zone_high.x = zones[i].zone_low.x + l0;
        zones[i].zone_high.y = zones[i].zone_low.y + l0;

        // printf("hindex: %d, low: %lf %lf, high: %lf %lf\n", i, zones[i].zone_low.x, zones[i].zone_low.y, zones[i].zone_high.x, zones[i].zone_high.y);
    }

    // printf("znx: %d, zny: %d\n", zone_num_x, zone_num_y);    

    // particle to zone
    for (int i = 0; i < num_particles; i++) {
        // printf("Particle: %d, posx: %lf, posy: %lf\n", i, bodies[i].px_, bodies[i].py_);

        int x = (bodies[i].px_ - pos_min.x) / l0;
        int y = (bodies[i].py_ - pos_min.y) / l0;

        // printf("Particle: %d, x: %d, y: %d\n", i, x, y);

        int hindex = d2h(x, y, zone_num_x, zone_num_y);
        // printf("hindex: %d push\n", hindex);
        zones[hindex].push_body(bodies[i]);
    }

    // printf("Particle to zone over.\n");    

    // print bodys
    // for (int i = 0; i < zone_num; i++) {
    //     for (int j = 0; j < z.size_; j++) {
    //         zones[i].bodys_[j].print();
    //     }
    // }

    // print zone size
    // for (int i = 0; i < zone_num_x * zone_num_y; i++) {
    //     printf("hindex: %d:", i);
    //     zones[i].print_zone_size();
    // }

    // only one core

    timeSt = dtime();
    for (int step = 0; step < cycles_times; step++) {
        // printf("-----------------Step: %d-----------------\n", step);
        for (int i = 0; i < zone_num; i++) {
            for (int part_index = 0; part_index < zones[i].size_; ++part_index) {
                zones[i].bodys_[part_index].ax_ = 0;
                zones[i].bodys_[part_index].ay_ = 0;
            }
            zones[i].calculate_a_self();

            // printf("My hindex: %d, self compute over.\n", i);

            // find neighbor
            // printf("znx: %d, zny: %d\n", zone_num_x, zone_num_y);   
            Pos<int> pos_t = h2d(i, zone_num_x, zone_num_y);
            // printf("My hindex: %d, dicard xy: %d %d\n", i, pos_t.x, pos_t.y);
            for (int n_j = 0; n_j < 8; n_j++) {
                int nei_x = pos_t.x + nei_dx[n_j];
                int nei_y = pos_t.y + nei_dy[n_j];

                if (nei_x < 0 || nei_x >= zone_num_x || nei_y < 0 || nei_y >= zone_num_y) continue;

                // printf("My hindex: %d, neighbor xy: %d %d, n_j: %d\n", i, nei_x, nei_y, n_j);
                int nei_hindex = d2h(nei_x, nei_y, zone_num_x, zone_num_y);

                // printf("My hindex: %d, will compute hindex: %d\n", i, nei_hindex);
                if (zones[nei_hindex].size_ > 0) {
                    // printf("In zone %d & zone %d\n", i, nei_hindex);
                    zones[i].calculate_a_neighbor(zones[nei_hindex]);
                }
            }
        }

        std::vector<Body> out_bodys;

        // update v & p
        for (int i = 0; i < zone_num; i++) {
            Zone* z = &zones[i];
            calculate_v_zone(z, max_v);
            calculate_p_zone(z, out_bodys);
        }

        // printf("Out bodys size: %d\n", out_bodys.size());

        // handle out bodys
        for (int i = 0; i < out_bodys.size(); i++) {
            if (out_bodys[i].px_ < pos_min.x || out_bodys[i].py_ < pos_min.y || 
                out_bodys[i].px_ > pos_max.x || out_bodys[i].py_ > pos_max.y) {
                // printf("Body out! Information:\n");
                // out_bodys[i].print();
            }
            else {
                int zone_x = (out_bodys[i].px_- pos_min.x) / l0;
                int zone_y = (out_bodys[i].py_- pos_min.y) / l0;

                int hindex = d2h(zone_x, zone_y, zone_num_x, zone_num_y);
                // printf("Zone: %d will into particle.\n", hindex);
                zones[hindex].push_body(out_bodys[i]);
            }
        }

        // // check point
        // if (step == 0/*cycles_times - 1*/)
        //     for (int i = 0; i < zone_num; i++) {
        //         Zone& z = zones[i];
        //         for (int j = 0; j < z.size_; j++) {
        //             z.bodys_[j].print();
        //         }
        //     }
    }

    timeSt = dtime() - timeSt;
    printf("All steps spend %.10lf without mpi.\n", timeSt);

    std::vector<Check_pos> check;
    for (int i = 0; i < zone_num; i++) {
        Zone& z = zones[i];
        for (int j = 0; j < z.size_; j++) {
            check.push_back(Check_pos(z.bodys_[j].px_, z.bodys_[j].py_));
        }
    }
    std::sort(check.begin(), check.end());
    FILE* outfile = fopen("mpi_no_result.txt", "w");
    if (outfile != NULL) {
        for (int i = 0; i < check.size(); i++) {
            fprintf(outfile, "%.4lf %.4lf\n", check[i].px, check[i].py);
        }
    }
    else {
        return 1;
    }

    fclose(file);

    MPI_Finalize();
    // athread_halt();

    return 0;
}

