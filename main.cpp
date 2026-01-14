//
// 2D Graphing by Interval Arithmetic, Implemented in C++
// References: https://zhuanlan.zhihu.com/p/1988701719459284501 and https://www.dgp.toronto.edu/public_user/mooncake/papers/SIGGRAPH2001_Tupper.pdf
// Copyright (c) 2026 YCX. Licensed under MIT License.
//

#include "interval.hpp"
#include <vector>
#include <utility>
#include <string>

int main() {
    Interval xrange(-10, 10);
    Interval yrange(-10, 10);

    int Nx = 2001;
    int Ny = 2001;

    CmpEqual cmp_eq;
    CmpLesser cmp_le;
    CmpGreater cmp_ge;

    // example 2D graphs
    std::vector<std::pair<Func2D, Cmp *> > curve_list = {
        // inexplicit curves
        { [](const IntervalSet & x, const IntervalSet & y) { return powOp(x, 3) + powOp(y, 3) - 15*x*y; }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(powOp(x, 2)) + sinOp(powOp(y, 2)) - 1.0; }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(y / x) + sinOp(x * y) - 1.0; }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return expOp(sinOp(x) + cosOp(y)) - sinOp(expOp(x + y)); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return (y - 5) * cosOp(4 * powOp(powOp(x - 4, 2) + powOp(y, 2), 0.5)) - x * sinOp(2 * powOp(powOp(x, 2) + powOp(y, 2), 0.5)); }, &cmp_ge },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(x * y) + sinOp(powOp(y, 2)) - 1.0; }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(4 * powOp(powOp(x, 2) + powOp(y, 2), 0.5)) - cosOp(x + y) + cosOp(x * powOp(y, 2)); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return powOp(absOp(x), 2.0/3.0) + powOp(absOp(y), 2.0/3.0) - 3; }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return x * y * sinOp(x * y) * cosOp(x * y) - (x + y); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(powOp(x, 2) + powOp(y, 2)) - cosOp(x * y); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(powOp(x, 2) + powOp(y, 2)) - cosOp(powOp(x, 2) - powOp(y, 2)); }, &cmp_le },
        { [](const IntervalSet & x, const IntervalSet & y) { return sinOp(powOp(x, 2) - powOp(y, 2)) - sinOp(x + y) - cosOp(x * y); }, &cmp_le },
        // explicit functions
        { [](const IntervalSet & x, const IntervalSet & y) { return y - sinOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - tanOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - 5 * sinOp(20 / x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - arccosOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - arctanOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - expOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - lnOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - absOp(x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - sqrtOp(x); }, &cmp_le },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - rootOp(x, -2); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - rootOp(x, 3); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - rootOp(x, -3); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - powOp(x, x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - powOp(x, 3); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - powOp(x, -2); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - powOp(x, -2.1); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - powOp(2, x); }, &cmp_eq },
        { [](const IntervalSet & x, const IntervalSet & y) { return y - 1.0 / (powOp(x, 2) - 3.0 * x + 2.0); }, &cmp_eq },
    };

    // plot & write file
    std::string out_dir = "./plot_data";
    size_t curve_id = 0;
    for (auto & p : curve_list) {
        Func2DPlotter plotter(p.first, p.second, xrange, yrange, Nx, Ny);
        plotter.plot();
        plotter.writeFile(out_dir + "/data_" + std::to_string(curve_id) + ".bin");
        curve_id++;
    }

    return 0;
}
