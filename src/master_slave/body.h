#ifndef _BODY_H_
#define _BODY_H_
#define real_type double

inline bool real_type_equal (real_type r1, real_type r2) {
    if (r1 - r2 <= 1e-5) return true;
    return false;
}

struct Packed_bodies {
    real_type px;
    real_type py;
    real_type m;
};

template <typename A>
struct Pos {
    A x;
    A y;
};

class Check_pos {
public: 
    real_type px;
    real_type py;
    Check_pos(real_type x, real_type y): px(x), py(y) {};
    bool operator<(const Check_pos& other) const {
        if (px == other.px) {
            return py < other.py;
        }
        else return px < other.px;
    }
};

struct BodySendAll {
    real_type px;
    real_type py;
    real_type vx;
    real_type vy;
    real_type m;
};

class Body {
public:
    Body(real_type px, real_type py, 
        real_type vx, real_type vy, real_type m = 1):
        px_(px), py_(py), vx_(vx), vy_(vy), m_(m) {};

    void print() {
        printf("px:%.4lf, py: %.4lf, vx: %.4lf, vy: %.4lf, ax: %.8lf, ay: %.8lf\n", px_, py_, vx_, vy_, ax_, ay_);
    }

    void print_result() {
        printf("px:%.4lf, py: %.4lf, vx: %.4lf, vy: %.4lf\n", px_, py_, vx_, vy_);
    }

    void print_result(int index) {
        printf("Zone: %d px:%.4lf, py: %.4lf, vx: %.4lf, vy: %.4lf\n", index, px_, py_, vx_, vy_);
    }

    bool operator== (Body& b1) {
        if (real_type_equal(b1.vx_, this->vx_) && real_type_equal(b1.vy_, this->vy_)) {
            if (real_type_equal(b1.px_, this->px_) && real_type_equal(b1.py_, this->py_)) {
                return true;
            }
        }
        else {
            printf("Not equal:\n");
            this->print();
            b1.print();
        }
        return false;
    }



// private:
    real_type px_;
    real_type py_;

    real_type vx_;
    real_type vy_;

    real_type ax_{0};
    real_type ay_{0};

    real_type m_;
};

#endif