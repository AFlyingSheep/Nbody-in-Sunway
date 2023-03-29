#ifndef _CALCULATE_H_
#define _CALCULATE_H_
#include <vector>
#include "zone.h"
#include "static.h"
#include "body.h"

class Zone;
class Body;

void calculate_v(Body&, real_type);
void calculate_v_zone(Zone*, real_type);
void calculate_v_zone(Zone*, size_t, size_t, real_type);
void calculate_p(Body&);
void calculate_p_zone(Zone*, std::vector<Body>&);

#endif