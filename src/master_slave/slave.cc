#include <crts.h>
#include <cmath>
#include "body.h"
#include "calculate.h"
// #include "main.h"
// #include "zone.h"

typedef struct {
    Body* body;
    int body_num;

} Calculate_a_zone_para;

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

extern "C" void calculate_a_cpe(void* cazp) {
    Calculate_a_zone_para* cazp_s = (Calculate_a_zone_para*)cazp;
    int body_num = cazp_s->body_num;
    Body* bodies = cazp_s->body;
    const int my_id = CRTS_tid;
    // // 测试一下是否接收到了？
    // printf("CPE ID: %d, body_num: %d, bodies: &%d.\n", my_id, body_num, bodies);

    // 根据自己的id判断计算bodies中的哪一部分
    int num_parts = 64;
    int part_size = (body_num + num_parts - 1) / num_parts;
    
    int start = my_id * part_size;
    if (start >= body_num) return;
    int end = (my_id + 1) * part_size;
    // printf("C num: %d.\n", end - start);
    // printf("Body num: %d\n", body_num);

    // printf("CPE ID: %d, start, end: %d, %d.\n", my_id, start, end);
    // 与其他body进行运算
    
    for (int now_part_index = start; now_part_index < end; now_part_index++) {
        bodies[now_part_index].ax_ = 0;
        bodies[now_part_index].ay_ = 0;
        for (int part_index = 0; part_index < body_num; part_index++) {
            if (now_part_index == part_index) continue;
            calculate_a(bodies[now_part_index], bodies[part_index]);
        }
    }
}

typedef struct {
    Body* body;
    Body* other_body;
    int body_num;
    int other_body_num;

} Calculate_a_zones_para;

extern "C" void calculate_a_zones_cpe(void* cazp) {
    Calculate_a_zones_para* cazp_s = (Calculate_a_zones_para*)cazp;
    int body_num = cazp_s->body_num;
    int other_body_num = cazp_s->other_body_num;
    Body* bodies = cazp_s->body;
    Body* other_bodies = cazp_s->other_body;

    const int my_id = CRTS_tid;
    // // 测试一下是否接收到了？
    // printf("CPE ID: %d, body_num: %d, bodies: &%d.\n", my_id, body_num, bodies);

    // 根据自己的id判断计算bodies中的哪一部分
    int num_parts = 64;
    int part_size = (body_num + num_parts - 1) / num_parts;
    
    int start = my_id * part_size;
    if (start >= body_num) return;
    int end = (my_id + 1) * part_size;
    
    for (int now_part_index = start; now_part_index < end; now_part_index++) {
        for (int part_index = 0; part_index < other_body_num; part_index++) {
            calculate_a(bodies[now_part_index], other_bodies[part_index]);
        }
    }
}

typedef struct {
    Body* body;
    Packed_bodies* other_body;
    int body_num;
    int other_body_num;

} Calculate_a_comm_zones_para;

extern "C" void calculate_a_comm_zones_cpe(void* cazp) {
    Calculate_a_comm_zones_para* cazp_s = (Calculate_a_comm_zones_para*)cazp;
    int body_num = cazp_s->body_num;
    int other_body_num = cazp_s->other_body_num;
    Body* bodies = cazp_s->body;
    Packed_bodies* other_bodies = cazp_s->other_body;

    const int my_id = CRTS_tid;
    // // 测试一下是否接收到了？
    // printf("CPE ID: %d, body_num: %d, bodies: &%d.\n", my_id, body_num, bodies);

    // 根据自己的id判断计算bodies中的哪一部分
    int num_parts = 64;
    int part_size = (body_num + num_parts - 1) / num_parts;
    
    int start = my_id * part_size;
    if (start >= body_num) return;
    int end = (my_id + 1) * part_size;
    
    for (int now_part_index = start; now_part_index < end; now_part_index++) {
        for (int part_index = 0; part_index < other_body_num; part_index++) {
            calculate_a(bodies[now_part_index], other_bodies[part_index]);
        }
    }
}