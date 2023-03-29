#include <crts.h>
#include <mpi.h>
#include <cmath>
// #include "main.h"
#include "zone.h"
#include "body.h"

extern "C" void SLAVE_FUN(calculate_a_cpe)(void*);
extern "C" void SLAVE_FUN(calculate_a_zones_cpe)(void*);
extern "C" void SLAVE_FUN(calculate_a_comm_zones_cpe)(void*);

// extern "C" int __real_athread_spawn(void*, void*, int);
// extern "C" int athread_join();
// extern "C" void calculate_a_cpe(void*);

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

void calculate_a(Body& ba, Packed_bodies& bb) {
    real_type dx = ba.px_ - bb.px, dy = ba.py_ - bb.py;
    // printf("dxdydz%lf%lf%lf", dx, dy, dz);
    real_type d_2 = dx * dx + dy * dy;
    real_type d = sqrt(d_2);

    // printf("%lf %lf %lf", d, d_2, dx);

    // printf("dxdydz%lf", d);
    // F = GMm / r3
    // Fx = GMm / r2 * (rx / r)
    // ax = GM / r2 * (rx / r)

    ba.ax_ -= G * bb.m / (pow(d, 3) + ITA) * dx;
    ba.ay_ -= G * bb.m / (pow(d, 3) + ITA) * dy;
}

typedef struct {
    Body* body;
    int body_num;

} Calculate_a_zone_para;

// test this
void calculate_a_zone(Zone* z) {
    int size = z->size_;
    // 将zone中的粒子分配给CPE进行计算
    const int cpe_num = 64;
    int per_body_of_CPE = std::floor(size / cpe_num);
    Calculate_a_zone_para cazp;
    cazp.body = z->bodys_;
    cazp.body_num = z->size_;

    // __real_athread_spawn(calculate_a_cpe, &cazp, 1);
    // printf("Bodies address:&%d", z->bodys_);
    athread_spawn((void*)(&SLAVE_FUN(calculate_a_cpe)), &cazp);
    athread_join();
    // for (size_t i = 0; i < z->size_; ++i) {
    //     z->bodys_[i].ax_ = 0;
    //     z->bodys_[i].ay_ = 0;
    //     // z->bodys_[i].az_ = 0;
    //     for (size_t j = 0; j < z->size_; ++j) {
    //         if (i == j) continue;
    //         calculate_a(z->bodys_[i], z->bodys_[j]);
    //     }
    // }
}

typedef struct {
    Body* body;
    Body* other_body;
    int body_num;
    int other_body_num;

} Calculate_a_zones_para;

void calculate_a_zones(Zone* z1, Zone* z2) {
    int z1_size = z1->size_;
    int z2_size = z2->size_;

    const int cpe_num = 64;
    int per_body_of_CPE = std::floor(z1_size / cpe_num);
    Calculate_a_zones_para cazp;
    cazp.body = z1->bodys_;
    cazp.body_num = z1->size_;
    cazp.other_body = z2->bodys_;
    cazp.other_body_num = z2->size_;

    athread_spawn((void*)(&SLAVE_FUN(calculate_a_zones_cpe)), &cazp);
    athread_join();
}

typedef struct {
    Body* body;
    Packed_bodies* other_body;
    int body_num;
    int other_body_num;

} Calculate_a_comm_zones_para;

void calculate_a_zones(Zone* z1, Comm_zone* z2) {
    int z1_size = z1->size_;
    int z2_size = z2->size_;

    const int cpe_num = 64;
    int per_body_of_CPE = std::floor(z1_size / cpe_num);
    Calculate_a_comm_zones_para cazp;
    cazp.body = z1->bodys_;
    cazp.body_num = z1->size_;
    cazp.other_body = z2->bodys_;
    cazp.other_body_num = z2->size_;

    athread_spawn((void*)(&SLAVE_FUN(calculate_a_comm_zones_cpe)), &cazp);
    athread_join();
}