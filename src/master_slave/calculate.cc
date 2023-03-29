#include <cmath>
#include <vector>
#include "calculate.h"

class Zone;
class Body;

void calculate_v(Body& b, real_type max_v) {
    b.vx_ += b.ax_ * delta_t;
    b.vy_ += b.ay_ * delta_t;
    // b.vz_ += b.az_ * delta_t;
    if (std::abs(b.vx_) > max_v) b.vx_ = (b.vx_ > 0)? max_v : -1.0 * max_v;
    if (std::abs(b.vy_) > max_v) b.vy_ = (b.vy_ > 0)? max_v : -1.0 * max_v;
}

void calculate_v_zone(Zone* z, real_type max_v) {
    for (size_t i = 0; i < z->size_; ++i) {
        calculate_v(z->bodys_[i], max_v);
    }
}

void calculate_v_zone(Zone* z, size_t start, size_t end, real_type max_v) {
    for (size_t i = start; i < end; ++i) {
        calculate_v(z->bodys_[i], max_v);
    }
}

void calculate_p(Body& b) {
    b.px_ += b.vx_ * delta_t;
    b.py_ += b.vy_ * delta_t;
    // b.pz_ += b.vz_ * delta_t;
}

void calculate_p_zone(Zone* z, std::vector<Body>& out_bodys) {
    for (size_t i = 0; i < z->size_; ++i) {
        Body& b = z->bodys_[i];
        b.px_ += b.vx_ * delta_t;
        b.py_ += b.vy_ * delta_t;

        // printf("Body px py: %lf, %lf.\n", b.px_, b.py_);
        // printf("zone_low_x, y, l0 = %lf %lf %lf.\n", z->zone_low.x, z->zone_low.y, z->l0);
        
        if (b.px_ < z->zone_low.x || b.py_ < z->zone_low.y || 
            b.px_ > z->zone_low.x + z->l0 || b.py_ > z->zone_low.y + z->l0) {
            out_bodys.push_back(b);
            z->pop(i);
        }
    }
}
