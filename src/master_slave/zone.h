#ifndef _ZONE_H_
#define _ZONE_H_
#include <stdio.h>
#include <string.h>
#define real_type double

#include "body.h"
#include "static.h"
#include "calculate.h"

class Comm_zone;

void calculate_a(Body&, Body&);
void calculate_a(Body&, Packed_bodies&);
void calculate_p(Body&);
void calculate_v(Body&);

class Zone {
public:
    Pos<real_type> zone_low;
    Pos<real_type> zone_high;

    real_type zone_height;
    real_type zone_width;

    Body* bodys_;
    real_type l0;
    size_t size_ {0};

    Zone() {
        bodys_ = (Body*)malloc(sizeof(Body) * 4);
        malloc_size = 4;
    }

    void print_zone_size() {
        printf("Zone size: %d\n", size_);
    }
    
    void push_body(const Body&);
    void calculate_a_self();
    void calculate_a_neighbor(Zone&);
    void calculate_a_neighbor(Comm_zone&);
    void pop(int);
    void print_body(int);

private:
    size_t malloc_size;
};

class Comm_zone {
public:
    Packed_bodies* bodys_;
    size_t size_ {0};

    Comm_zone() {
        bodys_ = (Packed_bodies*)malloc(sizeof(Packed_bodies) * 4);
        malloc_size = 4;
    }

    void print_zone_size() {
        printf("Zone size: %d\n", size_);
    }
    
    void push_body(const Packed_bodies&);
    void print_body();

private:
    size_t malloc_size;
};

#endif