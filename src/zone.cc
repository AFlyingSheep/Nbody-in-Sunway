#include <cstdio>
#include <iostream>
#include <string.h>

#include "zone.h"

void calculate_a(Body&, Body&);
void calculate_a(Body&, Packed_bodies&);
void calculate_p(Body&);
void calculate_v(Body&);

void Zone::push_body(const Body& body) {
    size_++;
    if (size_ > malloc_size) {
        malloc_size *= 2;
        bodys_ = (Body*)realloc(bodys_, malloc_size * sizeof(Body));
    }
    bodys_[size_ - 1] = body;
}

void Zone::pop(int index) {
    size_--;
    if (index == size_) return;
    bodys_[index] = bodys_[size_];
}

void Zone::calculate_a_self() {
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < size_; j++) {
            if (i == j) continue;
            calculate_a(bodys_[i], bodys_[j]);
        }
    }
}

void Zone::calculate_a_neighbor(Zone& z) {
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < z.size_; j++) {
            // printf("body: %d vs body: %d\n", i, j);
            calculate_a(bodys_[i], z.bodys_[j]);
        }
    }
}

void Zone::calculate_a_neighbor(Comm_zone& z) {
    for (int i = 0; i < size_; i++) {
        for (int j = 0; j < z.size_; j++) {
            calculate_a(bodys_[i], z.bodys_[j]);
        }
    }
}

void Zone::print_body(int index) {
    for (int i = 0; i < size_; i++) {
        bodys_[i].print_result(index);
    }
}

void Comm_zone::push_body(const Packed_bodies& packed_bodies) {
    size_++;
    if (size_ > malloc_size) {
        malloc_size *= 2;
        bodys_ = (Packed_bodies*)realloc(bodys_, malloc_size * sizeof(Body));
    }
    bodys_[size_ - 1] = packed_bodies;
}

void Comm_zone::print_body() {
    for (int i = 0; i < size_; i++) {
        printf("Packed bodies index: %d, px: %lf, py %lf, m %lf.\n", i, bodys_[i].px, bodys_[i].py, bodys_[i].m);
    }
}