#include <cstdio>
#include <math.h>
#include <cstdlib>
#include <mpi.h>
#include <assert.h>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <algorithm>
#include <crts.h>


#include "main.h"
#include "utils.h"
#include "zone.h"
#include "body.h"
#include "calculate.h"
#include "master.h"

#define FILE_NAME "particles2.txt"

// real_type delta_t = 1 / timestep;
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

struct Packed_bodies;

const int nei_dx[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
const int nei_dy[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

// extern void calculate_a(Body&, Body&);
// extern void calculate_a(Body&, Packed_bodies&);
extern void calculate_a_zone(Zone*);
// extern void calculate_a_zone(Zone*, size_t, size_t);
// extern void calculate_p(Body&);
// extern void calculate_v(Body&);
extern void calculate_p_zone(Zone*, std::vector<Body>&);
extern void calculate_v_zone(Zone*, real_type);

void calculate_test(const Zone* z, size_t start, size_t end, real_type dt) {
    for (size_t i = start; i < end; ++i) {
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
    for (int i = start; i < end; ++i) {
        z->bodys_[i].print_result();
    }
}

bool check(Zone* z1, Zone* z2) {
    // printf("z1: %d, z2: %d\n", z1->size_, z2->size_);
    assert(z1->size_ == z2->size_);

    // print_body(z2, true);
    for (size_t i = 0; i < z1->size_; ++i) {
        Body& b1 = z1->bodys_[i], &b2 = z2->bodys_[i];
        if (b1 == b2) continue;
        return false;
    }
    return true;
}

template <typename T>
struct Pos;
struct Packed_bodies;

int main(int argc, char** argv) {

    MPI_Init(NULL, NULL);
    int comm_sz;
    int my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    CRTS_init();

    printf("Rank: %d / %d says hello!\n", my_rank, comm_sz - 1);
    
    double timeSt;

    int zone_num_x, zone_num_y, zone_num;
    real_type l0;
    int* allocate_list = (int*)malloc(sizeof(int) * comm_sz);
    int* particle_num_in_zone;
    Zone* zones;
    Zone* l_zones;
    Pos<real_type> pos_max, pos_min;
    
    if (my_rank == 0) {
        // printf("-----------Data input----------\n");
        // input data
        FILE *file = fopen(FILE_NAME, "r");
        if (file == nullptr) {
            printf("Open file error!\n");
            return 1;
        }

        FILE *out_file = fopen("output.txt", "w");
        if (file == nullptr) {
            printf("Create file error!\n");
            return 1;
        }

        real_type x, y, vx, vy;
        int num_particles;

        fscanf(file, "%d %lf", &num_particles, &l0);
        Body* bodies = (Body*)malloc(sizeof(Body) * num_particles);

        pos_max.x = -1000;
        pos_max.y = -1000;
        pos_min.x = 1000;
        pos_min.y = 1000;

        for (int i = 0; i < num_particles; ++i) {
            fscanf(file, "%lf %lf %lf %lf", &x, &y, &vx, &vy);
            // printf("粒子 %d 的坐标为 (%lf, %lf)\n", i+1, x, y);
            bodies[i] = Body(x, y, vx, vy);
            // bodies[i].print_result();

            pos_max.x = std::max(pos_max.x, x);
            pos_max.y = std::max(pos_max.y, y);

            pos_min.x = std::min(pos_min.x, x);
            pos_min.y = std::min(pos_min.y, y);
        }

        fclose(file);
        fclose(out_file);

        pos_min.x -= l0 / 2;
        pos_min.y -= l0 / 2;
        pos_max.x += l0 / 2;
        pos_max.y += l0 / 2;

        real_type dx = pos_max.x - pos_min.x;
        real_type dy = pos_max.y - pos_min.y;

        // printf("dx: %lf, dy: %lf\n", dx, dy);

        // num of zones
        zone_num_x = std::ceil(dx / l0);
        zone_num_y = std::ceil(dy / l0);

        // **判断越界
        zone_num = zone_num_x * zone_num_y;
        // printf("Zone num: %d\n", zone_num);
        particle_num_in_zone = new int[zone_num] {0};

        // nums of every zone:
        // commsz有的没分配
        int per_zone_num = std::ceil(zone_num / comm_sz);
        // printf("Per zone num: %d\n", per_zone_num);

        // creat f list
        // just easy to allocate
        for (int i = 0; i < comm_sz; ++i) {
            allocate_list[i] = i * per_zone_num;
            // printf("Node: %d, list: %d\n", i, allocate_list[i]);
        }

        zones = (Zone*)malloc(sizeof(Zone) * zone_num_x * zone_num_y);

        for (int i = 0; i < zone_num; ++i) {
            zones[i] = Zone();
            Pos<int> p = h2d(i, zone_num_x, zone_num_y);
            zones[i].zone_low.x = p.x * l0 + pos_min.x;
            zones[i].zone_low.y = p.y * l0 + pos_min.y;

            zones[i].zone_high.x = zones[i].zone_low.x + l0;
            zones[i].zone_high.y = zones[i].zone_low.y + l0;

            // printf("hindex: %d, low: %lf %lf, high: %lf %lf\n", i, zones[i].zone_low.x, zones[i].zone_low.y, zones[i].zone_high.x, zones[i].zone_high.y);
        }

        printf("znx: %d, zny: %d\n", zone_num_x, zone_num_y);    

        // particle to zone
        for (int i = 0; i < num_particles; ++i) {
            // printf("Particle: %d, posx: %lf, posy: %lf\n", i, bodies[i].px_, bodies[i].py_);

            int x = (bodies[i].px_ - pos_min.x) / l0;
            int y = (bodies[i].py_ - pos_min.y) / l0;

            // printf("Particle: %d, x: %d, y: %d\n", i, x, y);

            int hindex = d2h(x, y, zone_num_x, zone_num_y);
            // printf("hindex: %d push\n", hindex);
            zones[hindex].push_body(bodies[i]);
            ++particle_num_in_zone[hindex];
        }

        // printf("Comm 0 Particle to zone over.\n");
        // for (int i = 0; i < zone_num; ++i) {
        //     if (zones[i].size_ > 0) printf("Zone: %d size is %d\n", i, zones[i].size_);
        // }

    }

    MPI_Bcast(allocate_list, comm_sz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zone_num_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zone_num_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zone_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&l0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (my_rank != 0) particle_num_in_zone = new int[zone_num] {0};
    MPI_Bcast(particle_num_in_zone, zone_num, MPI_INT, 0, MPI_COMM_WORLD);

    real_type pos_buf[2];
    pos_buf[0] = pos_min.x;
    pos_buf[1] = pos_min.y;

    MPI_Bcast(pos_buf, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    pos_min.x = pos_buf[0];
    pos_min.y = pos_buf[1];

    real_type max_v = l0 / delta_t;
    // printf("----!!!!posminxy: %lf, %lf.\n", pos_min.x, pos_min.y);

    int l_zone_num;
    if (my_rank != comm_sz - 1) l_zone_num = allocate_list[my_rank + 1] - allocate_list[my_rank];
    else l_zone_num = zone_num - allocate_list[my_rank];

    l_zones = new Zone[l_zone_num];

    // for (int i = 0; i < comm_sz; ++i) {
    //     printf("Node: %d, list: %d\n", i, allocate_list[i]);
    // }
    // printf("Node: %d, local zone num: %d\n", my_rank, l_zone_num);

    // print part_num_in_zone
    // if (my_rank == 0) {
    //     for (int i = 0; i < zone_num; ++i) {
    //         printf("%d, ", particle_num_in_zone[i]);
    //     }
    //     printf("\n");
    // }

    // send particles
    if (my_rank == 0) {
        for (int node_index = 1; node_index < comm_sz; ++node_index) {
            int send_particle_num = 0;
            int start = allocate_list[node_index];
            int end = (node_index == comm_sz - 1) ? zone_num : allocate_list[node_index + 1];
            // printf("Node index: %d, start: %d, end: %d.\n", node_index, start, end);
            for (int z_index = start; z_index < end; ++z_index) {
                send_particle_num += particle_num_in_zone[z_index];
            }
            real_type* sendbuf = new real_type[5 * send_particle_num];
            int sendbuf_p = 0;

            for (int z_index = start; z_index < end; ++z_index) {
                for (int p_index = 0; p_index < zones[z_index].size_; ++p_index) {
                    sendbuf[sendbuf_p++] = zones[z_index].bodys_[p_index].px_;
                    sendbuf[sendbuf_p++] = zones[z_index].bodys_[p_index].py_;
                    sendbuf[sendbuf_p++] = zones[z_index].bodys_[p_index].vx_;
                    sendbuf[sendbuf_p++] = zones[z_index].bodys_[p_index].vy_;
                    sendbuf[sendbuf_p++] = zones[z_index].bodys_[p_index].m_;
                }
            }
            assert(sendbuf_p == 5 * send_particle_num);
            // printf("I send 5 * %d to node %d.\n", send_particle_num, node_index);
            MPI_Send(sendbuf, 5 * send_particle_num, MPI_DOUBLE, node_index, 0, MPI_COMM_WORLD);
            delete[] sendbuf;
        }

        // itself
        for (int i = 0; i < l_zone_num; ++i) {
            Pos<int> p = h2d(i, zone_num_x, zone_num_y);
            l_zones[i].zone_low.x = p.x * l0 + pos_min.x;
            l_zones[i].zone_low.y = p.y * l0 + pos_min.y;
            l_zones[i].l0 = l0;

            for (int j = 0; j < particle_num_in_zone[i]; ++j) {
                l_zones[i].push_body(zones[i].bodys_[j]);
            }
        }
    }

    // recv particles
    else {
        // printf("----------recv test----------\n");
        int recv_particle_num = 0;
        int start = allocate_list[my_rank];
        int end = (my_rank == comm_sz - 1) ? zone_num : allocate_list[my_rank + 1];
        for (int z_index = start; z_index < end; ++z_index) {
            recv_particle_num += particle_num_in_zone[z_index];
        }
        MPI_Status status;
        int recv_count;
        int recv_point = 0;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &recv_count);
        assert(recv_count / 5 == recv_particle_num);
        
        real_type* recvbuf = new real_type[recv_count];
        MPI_Recv(recvbuf, recv_count, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int z_index = start; z_index < end; ++z_index) {
            int p_count = particle_num_in_zone[z_index];
            for (int p_index = 0; p_index < p_count; ++p_index) {
                l_zones[z_index - start].push_body(
                    Body(recvbuf[recv_point], recvbuf[recv_point + 1], 
                    recvbuf[recv_point + 2], recvbuf[recv_point + 3], 
                    recvbuf[recv_point + 4])
                );
                recv_point += 5;
            }
            // l_zones的左下角点&l0
            Pos<int> p = h2d(z_index, zone_num_x, zone_num_y);
            l_zones[z_index - start].zone_low.x = p.x * l0 + pos_min.x;
            l_zones[z_index - start].zone_low.y = p.y * l0 + pos_min.y;
            l_zones[z_index - start].l0 = l0;
        }
        // printf("recv_point: %d, recv_count: %d.\n", recv_point, recv_count);
        assert(recv_point == recv_count);
    }

    // print bodies
    int particle_size = 0;
    for (int i = allocate_list[my_rank]; i < ((my_rank == comm_sz - 1) ? zone_num : allocate_list[my_rank + 1]); i++) {
        // l_zones[i - allocate_list[my_rank]].print_body(i);
        particle_size += l_zones[i - allocate_list[my_rank]].size_;
        // printf("-------------Zone: %d, size: %d.\n", i, l_zones[i - allocate_list[my_rank]].size_);
    }
    printf("Node: %d, I have %d bodies.\n", my_rank, particle_size);
        
    int start = allocate_list[my_rank];
    int end = (my_rank == comm_sz - 1) ? zone_num : allocate_list[my_rank + 1];

    // printf("--Node index: %d, start: %d, end: %d.\n", my_rank, start, end);
    MPI_Barrier(MPI_COMM_WORLD);

    timeSt = dtime();
    double a_time = 0, pv_time = 0;
    double timeSt2;

    double times[10] = {0};


    for (int step = 0; step < cycles_times; ++step) {
        // printf("-----------------Step: %d-----------------\n", step);
        std::map<int, std::set<int>> send_list;
        std::map<int, int> my_neighbor_node;

        double temp = dtime();

        // 更新zone内部相互作用粒子加速度
        for (int zone_index = start; zone_index < end; ++zone_index) {
            // printf("Zone index: %d, start to calculate itself.\n", zone_index);
            timeSt2 = dtime();
            // double ta = dtime();

            // l_zones[zone_index - start].calculate_a_self();
            // ta = dtime() - ta;
            // if (step == 0) printf("Time no cpe: %lf.\n", ta);
            // ta = dtime();

            calculate_a_zone(&l_zones[zone_index - start]);
            // ta = dtime() - ta;
            // if (step == 0) printf("Time cpe: %lf.\n", ta);
            a_time += dtime() - timeSt2;
            // find neihber
            Pos<int> pos_t = h2d(zone_index, zone_num_x, zone_num_y);
            // printf("Zone index: %d, my pos: %d, %d.\n", zone_index, pos_t.x, pos_t.y);

            for (int n_j = 0; n_j < 8; n_j++) {
                int nei_x = pos_t.x + nei_dx[n_j];
                int nei_y = pos_t.y + nei_dy[n_j];

                if (nei_x < 0 || nei_x >= zone_num_x || nei_y < 0 || nei_y >= zone_num_y) continue;

                int nei_hindex = d2h(nei_x, nei_y, zone_num_x, zone_num_y);
                int nei_node = find_node_by_hindex(nei_hindex, allocate_list, comm_sz);

                // printf("Zone index: %d, neighbor hindex: %d in node %d.\n", zone_index, nei_hindex, nei_node);

                // 发送接受思路：
                // 发送list数据结构: map<int, set<int>> --> map<node index, set<zone index>>
                // zone_index代表欲发送到隔壁node的zone
                
                // 对当前的zone寻找邻居zone，看一看邻居zone位于哪一个node，
                // 如果在其他node(other_node_index)，则list添加: map[other_node_index].insert(this_zone_index);

                if (nei_node == my_rank) {
                    // 当zone刚好处在当前节点内，直到计算的时候拿过来就行了
                }
                else {
                    send_list[nei_node].insert(zone_index);
                    if (my_neighbor_node.find(nei_node) == my_neighbor_node.end()) {
                        my_neighbor_node.insert(std::make_pair(nei_node, 0));
                    }
                }
            }   // neighbor end
        }   // zones in node end.

        // // check send list
        // auto iter = send_list.begin();
        // for (; iter != send_list.end(); iter++) {
        //     printf("Node: %d <--", iter->first);
        //     for (auto iter2 = iter->second.begin(); iter2 != iter->second.end(); iter2++) {
        //         printf("%d ", *iter2);
        //     }
        //     printf("\n");
        // }

        MPI_Barrier(MPI_COMM_WORLD);
        times[0] += dtime() - temp;
        temp = dtime();

        // 根据sendlist来向指定节点发送数据
        for (auto iter = send_list.begin(); iter != send_list.end(); iter++) {
            int des_node_index = iter->first;
            std::set<int>& bodies_index = iter->second;

            // 根据send_list发送数据都来自哪个zone的
            int* send_index_list = new int[bodies_index.size()];
            // 每个zone含有多少个粒子
            int* send_count_list = new int[bodies_index.size()];

            // 包装需要发送的数据
            int send_length = 0;
            int count_point = 0;
            for (auto set_iter = bodies_index.begin(); set_iter != bodies_index.end(); ++set_iter) {
                send_length += l_zones[*set_iter - start].size_;
                send_count_list[count_point++] = l_zones[*set_iter - start].size_;
            }
            send_length *= 3;   // px py m

            // printf("Node: %d, send length: %d.\n", my_rank, send_length);

            double* send_buf = new double[send_length];
            int point = 0;
            int index = 0;
            for (auto set_iter = bodies_index.begin(); set_iter != bodies_index.end(); ++set_iter) {
                // 遍历每个粒子并包装
                send_index_list[index++] = *set_iter;
                for (int i = 0; i < l_zones[*set_iter - start].size_; i++) {
                    Body* bodies = l_zones[*set_iter - start].bodys_;
                    send_buf[point++] = bodies[i].px_;
                    send_buf[point++] = bodies[i].py_;
                    send_buf[point++] = bodies[i].m_;
                }
            }

            assert(point == send_length);
            assert(index == bodies_index.size());
            assert(count_point == index);

            // 异步发送
            MPI_Request request;
            MPI_Isend(send_index_list, index, MPI_INT, des_node_index,
                        0, MPI_COMM_WORLD, &request);
            MPI_Isend(send_buf, send_length, MPI_DOUBLE, des_node_index,
                        1, MPI_COMM_WORLD, &request);
            MPI_Isend(send_count_list, count_point, MPI_INT, des_node_index,
                        2, MPI_COMM_WORLD, &request);

            // TAG == 0 index_list
            // TAG == 1 bodies_data
        }

        // 转为接受，接受邻居发送过来的数据
        // 发送给多少个node，就需要接受多少个node的值
        std::set<int> recv_zone_index;
        int recv_node_count = send_list.size();
        Zone* recv_z;

        int err;
        MPI_Status status, status_next;
        std::vector<real_type> recv_buf;
        std::vector<int> recv_index_list;
        std::vector<int> recv_count_list;

        while (recv_node_count > 0) {

            std::vector<real_type> recv_buf_temp;
            std::vector<int> recv_index_list_temp;
            std::vector<int> recv_count_list_temp;

            int count;
            // 先测试接受send_index_list
            err = MPI_Probe(MPI_ANY_SOURCE, 0,
                            MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &count);
            recv_index_list_temp.resize(count);
            MPI_Recv(recv_index_list_temp.data(), count, MPI_INT, status.MPI_SOURCE, 
                        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("Node: %d, I recive list(data size: %d) from node %d.\n", my_rank, count, status.MPI_SOURCE);

            err = MPI_Probe(status.MPI_SOURCE, 1,
                            MPI_COMM_WORLD, &status_next);
            MPI_Get_count(&status_next, MPI_DOUBLE, &count);
            recv_buf_temp.resize(count);
            MPI_Recv(recv_buf_temp.data(), count, MPI_DOUBLE, status.MPI_SOURCE, 
                        1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // printf("Node: %d, I recive bodies(data size: %d) from node %d.\n", my_rank, count, status.MPI_SOURCE);

            err = MPI_Probe(status.MPI_SOURCE, 2,
                            MPI_COMM_WORLD, &status_next);
            MPI_Get_count(&status_next, MPI_INT, &count);
            recv_count_list_temp.resize(count);
            MPI_Recv(recv_count_list_temp.data(), count, MPI_INT, status.MPI_SOURCE, 
                        2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // test recv data
            assert(recv_count_list_temp.size() == recv_index_list_temp.size());
            for (int i = 0; i < recv_count_list_temp.size(); i++) {
                // printf("Recv particles from zone %d and size is %d.\n", recv_index_list_temp[i], recv_count_list_temp[i]);
            }
            // assert have right particle data
            int recv_length = 0;
            for (int i = 0; i < recv_count_list_temp.size(); i++) {
                recv_length += recv_count_list_temp[i];
            }
            assert(recv_length * 3 == recv_buf_temp.size());

            // 此处可以用reserve进一步提高性能
            recv_buf.insert(recv_buf.end(), recv_buf_temp.begin(), recv_buf_temp.end());
            recv_index_list.insert(recv_index_list.end(), recv_index_list_temp.begin(), recv_index_list_temp.end());
            recv_count_list.insert(recv_count_list.end(), recv_count_list_temp.begin(), recv_count_list_temp.end());

            --recv_node_count;
        }

        // 接收完毕，recv_index_list记录了zone号，recv_buf是接受的数据，recv_count_list记录对应zone的粒子数
        // 计算加速度
        // if (my_rank == 0) printf("----------Calculate with neighbor----------\n");
        // if (my_rank == 2) {
        //     for (int i = 0; i < recv_index_list.size(); i++) {
        //         printf("%d ", recv_index_list[i]);
        //     }
        //     printf("\n");
        // }
        // if (my_rank == 2) {
        //     for (int i = 0; i < recv_count_list.size(); i++) {
        //         printf("%d ", recv_count_list[i]);
        //     }
        //     printf("\n");
        // }
        // if (my_rank == 2) {
        //     for (int i = 0; i < recv_buf.size(); i++) {
        //         printf("%lf ", recv_buf[i]);
        //     }
        //     printf("\n");
        // }
        // if (my_rank == 2)

        MPI_Barrier(MPI_COMM_WORLD);
        times[1] += dtime() - temp;
        temp = dtime();

        for (int zone_index = start; zone_index < end; ++zone_index) {
            // find neihber
            Pos<int> pos_t = h2d(zone_index, zone_num_x, zone_num_y);
            // printf("Zone index: %d, my pos: %d, %d.\n", zone_index, pos_t.x, pos_t.y);
            for (int n_j = 0; n_j < 8; n_j++) {
                int nei_x = pos_t.x + nei_dx[n_j];
                int nei_y = pos_t.y + nei_dy[n_j];

                if (nei_x < 0 || nei_x >= zone_num_x || nei_y < 0 || nei_y >= zone_num_y) continue;

                int nei_hindex = d2h(nei_x, nei_y, zone_num_x, zone_num_y);
                // printf("My neighbor hindex is %d.\n", nei_hindex);
                int nei_node = find_node_by_hindex(nei_hindex, allocate_list, comm_sz);
                timeSt2 = dtime();
                if (nei_node == my_rank) {
                    l_zones[zone_index - start].calculate_a_neighbor(l_zones[nei_hindex - start]);
                    // calculate_a_zones(&l_zones[zone_index - start], &l_zones[nei_hindex - start]);
                }
                else {
                    // 寻找buf中所需zone开始的位置
                    int index_start = 0;
                    int find_index = 0;
                    for (; find_index < recv_index_list.size(); ++find_index) {
                        if (recv_index_list[find_index] == nei_hindex) break;
                        index_start += recv_count_list[find_index] * 3;
                    }
                    assert(find_index != recv_index_list.size());
                    Comm_zone temp_zone;
                    // printf("Start index: %d.\n", index_start);
                    for (int body_index = 0; body_index < recv_count_list[find_index]; ++body_index) {
                        Packed_bodies temp_body;
                        temp_body.px = recv_buf[index_start + body_index * 3];
                        temp_body.py = recv_buf[index_start + body_index * 3 + 1];
                        temp_body.m = recv_buf[index_start + body_index * 3 + 2];

                        temp_zone.push_body(temp_body);
                    }

                    l_zones[zone_index - start].calculate_a_neighbor(temp_zone);
                    // calculate_a_zones(&l_zones[zone_index - start], &temp_zone);
                }
                a_time += dtime() - timeSt2;
            }   // neighbor end
        }   // zones in node end.
        MPI_Barrier(MPI_COMM_WORLD);

        // printf("Node %d start to calculate v & p.\n", my_rank);

        times[2] += dtime() - temp;
        temp = dtime();

        // 计算速度和位置，如果遇到溢出的粒子则进行处理
        std::vector<Body> out_bodys;
        timeSt2 = dtime();
        for (int zone_index = start; zone_index < end; ++zone_index) {
            calculate_v_zone(&l_zones[zone_index - start], max_v);
            calculate_p_zone(&l_zones[zone_index - start], out_bodys);
        }
        pv_time += dtime() - timeSt2;
        for (Body& body: out_bodys) {
            if (body.px_ < pos_min.x || body.py_ < pos_min.y || 
                body.px_ > pos_min.x + zone_num_x * l0 || body.py_ > pos_min.y + zone_num_y * l0) {
                // printf("Body out! Information:\n");
                // body.print();
            }
            else {
                int zone_x = (body.px_ - pos_min.x) / l0;
                int zone_y = (body.py_ - pos_min.y) / l0;
                int body_hindex = d2h(zone_x, zone_y, zone_num_x, zone_num_y);
                int body_node = find_node_by_hindex(body_hindex, allocate_list, comm_sz);

                if (body_node == my_rank) {
                    l_zones[body_hindex - start].push_body(body);
                }
                else {
                    // 发送粒子
                    // 这里暂时先一个一个发送，性能可以提升
                    // 断言：溢出的粒子一定到了邻居的node
                    if (my_neighbor_node.find(body_node) == my_neighbor_node.end()) {
                        printf("Node: %d, out body node: %d.\n", my_rank, body_node);
                        printf("Out body's v = %lf %lf.\n", body.vx_, body.vy_);
                        assert(my_neighbor_node.find(body_node) != my_neighbor_node.end());
                    }
                    
                    ++my_neighbor_node[body_node];

                    // 包装粒子
                    double* packed_body = new double[6] {body.px_, body.py_, body.vx_, body.vy_, body.m_, body_hindex * 1.0};    
                    // px py vx vy m dest_zone

                    // 发送粒子
                    MPI_Request request;
                    // printf("Node: %d send a body, %lf %lf %lf %lf %lf to zone %d in node %d.\n", my_rank, 
                    //         packed_body[0], packed_body[1], packed_body[2], packed_body[3], packed_body[4], 
                    //         body_hindex, body_node);
                    MPI_Isend(packed_body, 6, MPI_DOUBLE, body_node, 3, MPI_COMM_WORLD, &request);
                }
            }
        }
        // 先把my_neighbor_node发送出去，确定需要接受多少粒子
        // check my_neighbor_node
        for (auto it : my_neighbor_node) {
            // printf("Node: %d has %d bodies send to node %d.\n", my_rank, it.second, it.first);
        }

        for (auto it : my_neighbor_node) {
            MPI_Request request;
            MPI_Isend(&it.second, 1, MPI_INT, it.first, 4, MPI_COMM_WORLD, &request);
        }
        int recv_out_body_count = 0;
        int recv_index_count = 0;
        for (auto it : my_neighbor_node) {
            int temp = 0;
            MPI_Recv(&temp, 1, MPI_INT, it.first, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            recv_out_body_count += temp;
        }
        // printf("Node: %d, recv count = %d.\n", my_rank, recv_out_body_count);

        // 然后接受粒子直到达到需要的数量
        while (recv_out_body_count > 0) {
            std::unique_ptr<real_type[]> temp_ptr(new real_type[6]);
            MPI_Status status;
            int count;
            MPI_Probe(MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &count);
            assert(count == 6);
            MPI_Recv(temp_ptr.get(), 6, MPI_DOUBLE, status.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // printf("Node: %d recv a body, %lf %lf %lf %lf %lf to zone %.0lf.\n", my_rank, 
            //                 temp_ptr[0], temp_ptr[1], temp_ptr[2], temp_ptr[3], temp_ptr[4], temp_ptr[5]);

            // 插入到zone中
            int insert_zone_index = static_cast<int>(temp_ptr[5]);
            l_zones[insert_zone_index - start].push_body(
                Body(temp_ptr[0], temp_ptr[1], temp_ptr[2], temp_ptr[3], temp_ptr[4])
            );
            // printf("Zone %d (start: %d) in node %d has insert a body.\n", insert_zone_index, start, my_rank);
            --recv_out_body_count;
        }
        // 迁移完成

        MPI_Barrier(MPI_COMM_WORLD);

        // 清除加速度
        for (int zone_index = start; zone_index < end; ++zone_index) {
            Zone& z = l_zones[zone_index - start];
            for (int body_index = 0; body_index < z.size_; ++body_index) {
                z.bodys_[body_index].ax_ = 0;
                z.bodys_[body_index].ay_ = 0;
            }
        }

        times[3] += dtime() - temp;
        temp = dtime();
    }

    // TAG3: 向邻居发送粒子
    // TAG4: 向邻居发送需要接受的粒子数量
    // TAG5: CHECKPOINT
    MPI_Barrier(MPI_COMM_WORLD);
    timeSt = dtime() - timeSt;
    if (my_rank == 0) printf("All steps spend %.10lf with mpi.\n", timeSt);
    if (my_rank == 0) printf("A steps spend %.10lf with mpi.\n", a_time);
    if (my_rank == 0) printf("P&V steps spend %.10lf with mpi.\n", pv_time);

    if (my_rank == 0) {
        for (int i = 0; i < 4; i++) printf("Steps %d spend %.10lf.\n", i, times[i]);
    }
    // check point
    std::vector<real_type> send_point;
    for (int zone_index = start; zone_index < end; zone_index++) {
        for (int part_index = 0; part_index < l_zones[zone_index - start].size_; part_index++) {
            send_point.push_back(l_zones[zone_index - start].bodys_[part_index].px_);
            send_point.push_back(l_zones[zone_index - start].bodys_[part_index].py_);
        }
    }

    if (my_rank == 0) {
        for (int node_index = 1; node_index < comm_sz; ++node_index) {
            MPI_Status status;
            int count;
            MPI_Probe(MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &count);
            std::vector<real_type> recv_point(count);

            MPI_Recv(recv_point.data(), count, MPI_DOUBLE, status.MPI_SOURCE, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            send_point.insert(send_point.end(), recv_point.begin(), recv_point.end());
        }
    }
    else {
        MPI_Send(send_point.data(), send_point.size(), MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
    }

    if (my_rank == 0) {
        std::vector<Check_pos> check;
        for (int i = 0; i < send_point.size() / 2; i++) {
            check.push_back(Check_pos(send_point[i * 2], send_point[i * 2 + 1]));
        }
        std::sort(check.begin(), check.end());
        FILE* outfile = fopen("mpi_result.txt", "w");
        if (outfile != NULL) {
            for (int i = 0; i < check.size(); i++) {
                fprintf(outfile, "%.4lf %.4lf\n", check[i].px, check[i].py);
            }
        }
        else {
            return 1;
        }
    }

    free(allocate_list);
    MPI_Finalize();
    athread_halt();

    return 0;
}

// int main(int argc, char** argv) {

//     MPI_Init(NULL, NULL);
//     int comm_sz;
//     int my_rank;
//     MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
//     MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

//     CRTS_init();

//     printf("Rank: %d / %d says hello!\n", my_rank, comm_sz - 1);
    
//     Zone z;
//     z.push_body(Body(0, 0, 0, 0, 1));
//     z.push_body(Body(0, 1, 0, 0, 1));
//     z.push_body(Body(1, 0, 0, 0, 1));
//     z.push_body(Body(1, 1, 0, 0, 1));

//     calculate_a_zone(&z);

//     MPI_Finalize();
//     athread_halt();

//     return 0;
// }