#include "../interval.hpp"

#include <cassert>
#include <iostream>

int main() {
    CmpEqual cmp_eq;
    Func2D logarithm = [](const IntervalSet & x, const IntervalSet & y) {
        return y - lnOp(x);
    };

    Func2DPlotter no_subpixels(logarithm, &cmp_eq, Interval(-10.0, 10.0), Interval(-10.0, 10.0), 201, 201);
    no_subpixels.plot(0);

    Func2DPlotter plotter(logarithm, &cmp_eq, Interval(-10.0, 10.0), Interval(-10.0, 10.0), 201, 201);
    plotter.plot(16);

    // The x=0 boundary pixel covers positive x values down to zero.  Therefore
    // it contains ln(x) for every visible negative y value, including values
    // far below ln(half a pixel), where the old implementation left a gap.
    int boundary_x = 100;
    int y_minus_ten = 0;
    int y_minus_six = 40;
    assert(no_subpixels.data[y_minus_ten + boundary_x * no_subpixels.Ny] == 0);
    assert(plotter.data[y_minus_ten + boundary_x * plotter.Ny] == 1);
    assert(plotter.data[y_minus_six + boundary_x * plotter.Ny] == 1);

    // A pixel wholly on the negative side remains outside ln(x)'s domain.
    int negative_x = 99;
    assert(plotter.data[y_minus_ten + negative_x * plotter.Ny] == 0);

    std::cout << "plotter regression tests passed\n";
    return 0;
}
