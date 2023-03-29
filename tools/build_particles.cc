#include <iostream>
#include <fstream>
#include <random>
using namespace std;

int main() {
    // 输入边界和粒子数量
    double x_min, x_max, y_min, y_max;
    int num_particles;
    double l0;
    cout << "x_min x_max y_min y_max：\n";
    cin >> x_min >> x_max >> y_min >> y_max;
    cout << "num of particles\n";
    cin >> num_particles;

    cout << "l0:\n";

    cin >> l0;

    // 随机生成粒子坐标
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis_x(x_min, x_max);
    uniform_real_distribution<double> dis_y(y_min, y_max);
    vector<pair<double, double>> particles;
    for (int i = 0; i < num_particles; i++) {
        double x = dis_x(gen);
        double y = dis_y(gen);
        particles.push_back(make_pair(x, y));
    }

    // 输出到文件
    ofstream outfile("particles.txt");
    outfile << num_particles << " " << l0 << endl;
    for (int i = 0; i < num_particles; i++) {
        outfile << particles[i].first << " " << particles[i].second << " 0 0"<< endl;
    }
    outfile.close();

    cout << "已生成 " << num_particles << " 个粒子的坐标，并保存到 particles.txt 中。" << endl;

    return 0;
}