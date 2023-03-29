#ifndef _MAIN_H_
#define _MAIN_H_

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#define real_type double
#include "body.h"
#include "static.h"

// const real_type G = 6.67e-5;

// const real_type timestep = 0.1;
// const int cycles_times = 2;
// const int ITA = 0.001;

class Body;
struct Packed_bodies;
template<typename A> struct Pos;
class Check_pos;
struct BodySendAll;

void calculate_a(Body&, Body&);
void calculate_a(Body&, Packed_bodies&);
void calculate_p(Body&);
void calculate_v(Body&);

int d2h(int x, int y, int zone_num_x, int zone_num_y) {
    int ans;
    bool flag = (y % 2 == 0);
    if (flag)   ans = zone_num_x * y + x;
    else        ans = zone_num_x * (y + 1) - x - 1;

    return ans;
}

Pos<int> h2d(int h, int zone_num_x, int zone_num_y) {
    Pos<int> ans;
    int y_index = std::floor(h / zone_num_x);
    if (y_index % 2 == 0) {
        ans.x = h % zone_num_x;
        ans.y = y_index;
    }
    else {
        ans.x = zone_num_x - h % zone_num_x - 1;
        ans.y = y_index;
    }

    return ans;
}

int find_node_by_hindex(int hindex, int* allocate_list, int size) {
    for (int i = 0; i < size - 1; i++) {
        if (hindex < allocate_list[i + 1]) return i;
    }
    return size - 1;
}

#endif
