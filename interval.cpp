//
// 2D Graphing by Interval Arithmetic, Implemented in C++
// References: https://zhuanlan.zhihu.com/p/1988701719459284501 and https://www.dgp.toronto.edu/public_user/mooncake/papers/SIGGRAPH2001_Tupper.pdf
// Copyright (c) 2026 YCX. Licensed under MIT License.
//

#include "interval.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

// addition, subtraction, multiplication, division
IntervalSet addOpInterval(const Interval & a, const Interval & b) {
    return Interval(a.lower + b.lower, a.upper + b.upper, a.well_defined && b.well_defined);
}
IntervalSet subOpInterval(const Interval & a, const Interval & b) {
    return Interval(a.lower - b.upper, a.upper - b.lower, a.well_defined && b.well_defined);
}
IntervalSet negOpInterval(const Interval & a) {
    return Interval(-a.upper, -a.lower, a.well_defined);
}
IntervalSet mulOpInterval(const Interval & a, const Interval & b) {
    return Interval(
        std::min({a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper}),
        std::max({a.lower * b.lower, a.lower * b.upper, a.upper * b.lower, a.upper * b.upper}),
        a.well_defined && b.well_defined
    );
}
IntervalSet divOpInterval(const Interval & a, const Interval & b) {
    IntervalSet b_inv;
    if (b.lower > 0.0 || b.upper < 0.0) {
        b_inv = Interval(1.0 / b.upper, 1.0 / b.lower, b.well_defined);
    } else if (b.lower < 0.0 && b.upper > 0.0) {
        b_inv = IntervalSet({
            Interval::toNegInf(1.0 / b.lower, false),
            Interval::toPosInf(1.0 / b.upper, false),
        });
    } else if (b.lower == 0.0) {
        b_inv = Interval::toPosInf(1.0 / b.upper, false);
    } else if (b.upper == 0.0) {
        b_inv = Interval::toNegInf(1.0 / b.lower, false);
    }
    return b_inv * a;
}

bool IntervalSet::isWellDefined() const {
    if (this->isEmpty()) {
        return false;
    }
    bool well_defined = true;
    for (const auto & intv : this->intervals) {
        well_defined = well_defined && intv.well_defined;
    }
    return well_defined;
}

double IntervalSet::getLowerBound() const {
    return std::min_element(
        this->intervals.begin(), this->intervals.end(),
        [](const Interval & a, const Interval & b) { return a.lower < b.lower; }
    )->lower;
}

double IntervalSet::getUpperBound() const {
    return std::max_element(
        this->intervals.begin(), this->intervals.end(),
        [](const Interval & a, const Interval & b) { return a.upper < b.upper; }
    )->upper;
}

// merge intervals that overlap
IntervalSet IntervalSet::mergeOverlapping() const {
    if (intervals.empty()) return *this;
    
    std::vector<Interval> sorted = intervals;
    std::sort(
        sorted.begin(), sorted.end(), 
        [](const Interval& a, const Interval& b) { return a.lower < b.lower; }
    );
    
    std::vector<Interval> merged;
    Interval current = sorted[0];
    
    for (size_t i = 1; i < sorted.size(); ++i) {
        if (current.upper >= sorted[i].lower) {
            current.upper = std::max(current.upper, sorted[i].upper);
        } else {
            merged.push_back(current);
            current = sorted[i];
        }
    }
    merged.push_back(current);
    
    return IntervalSet(merged);
}

IntervalSet IntervalSet::applyUnaryOp(UnaryOp func) const {
    std::vector<Interval> result;
    for (const auto & interval : intervals) {
        IntervalSet transformed = func(interval);
        auto transformedIntervals = transformed.getIntervals();
        result.insert(result.end(), transformedIntervals.begin(), transformedIntervals.end());
    }
    return IntervalSet(result);
}

IntervalSet IntervalSet::applyBinaryOp(const IntervalSet & other, BinaryOp func) const {
    std::vector<Interval> result;
    for (const auto & interval1 : intervals) {
        for (const auto & interval2 : other.getIntervals()) {
            IntervalSet transformed = func(interval1, interval2);
            auto transformedIntervals = transformed.getIntervals();
            result.insert(result.end(), transformedIntervals.begin(), transformedIntervals.end());
        }
    }
    return IntervalSet(result);
}

// exponent, logrithm, (arc) triangular functions, absolute value, n-th root, power
IntervalSet expOpInterval(const Interval & interval) {
    return Interval(exp(interval.lower), exp(interval.upper), interval.well_defined);
}
IntervalSet lnOpInterval(const Interval & interval) {
    if (interval.lower > 0.0) {
        return Interval(log(interval.lower), log(interval.upper), interval.well_defined);
    } else if (interval.upper < 0.0) {
        return IntervalSet();
    } else {
        return Interval::toNegInf(log(interval.upper), false);
    }
}
IntervalSet sinOpInterval(const Interval & interval) {
    double v1 = sin(interval.lower);
    double v2 = sin(interval.upper);
    double vmax = std::max(v1, v2);
    double vmin = std::min(v1, v2);
    int il = ceil(interval.lower / M_PI_2);
    int iu = floor(interval.upper / M_PI_2);
    int diff = iu - il;
    il = il % 4;
    if (il < 0) {
        il += 4;
    }
    iu = il + diff;
    if (il <= 1 && 1 <= iu) {
        vmax = 1.0;
    }
    if (il <= 3 && 3 <= iu) {
        vmin = -1.0;
    }
    if (il <= 5 && 5 <= iu) {
        vmax = 1.0;
    }
    if (il <= 7 && 7 <= iu) {
        vmin = -1.0;
    }
    return Interval(vmin, vmax, interval.well_defined);
}
IntervalSet cosOpInterval(const Interval & interval) {
    double v1 = cos(interval.lower);
    double v2 = cos(interval.upper);
    double vmax = std::max(v1, v2);
    double vmin = std::min(v1, v2);
    int il = ceil(interval.lower / M_PI_2);
    int iu = floor(interval.upper / M_PI_2);
    int diff = iu - il;
    il = il % 4;
    if (il < 0) {
        il += 4;
    }
    iu = il + diff;
    if (il <= 0 && 0 <= iu) {
        vmax = 1.0;
    }
    if (il <= 2 && 2 <= iu) {
        vmin = -1.0;
    }
    if (il <= 4 && 4 <= iu) {
        vmax = 1.0;
    }
    if (il <= 6 && 6 <= iu) {
        vmin = -1.0;
    }
    return Interval(vmin, vmax, interval.well_defined);
}
IntervalSet arcsinOpInterval(const Interval & interval) {
    if (interval.upper < -1.0 || interval.lower > 1.0) {
        return IntervalSet();
    }
    Interval interval0;
    interval0.well_defined = interval.well_defined;
    if (interval.lower < -1.0) {
        interval0.lower = -1.0;
        interval0.well_defined = false;
    } else {
        interval0.lower = interval.lower;
    }
    if (interval.upper > 1.0) {
        interval0.upper = 1.0;
        interval0.well_defined = false;
    } else {
        interval0.upper = interval.upper;
    }
    return Interval(asin(interval0.lower), asin(interval0.upper), interval0.well_defined);
}
IntervalSet arccosOpInterval(const Interval & interval) {
    if (interval.upper < -1.0 || interval.lower > 1.0) {
        return IntervalSet();
    }
    Interval interval0;
    interval0.well_defined = interval.well_defined;
    if (interval.lower < -1.0) {
        interval0.lower = -1.0;
        interval0.well_defined = false;
    } else {
        interval0.lower = interval.lower;
    }
    if (interval.upper > 1.0) {
        interval0.upper = 1.0;
        interval0.well_defined = false;
    } else {
        interval0.upper = interval.upper;
    }
    return Interval(acos(interval0.upper), acos(interval0.lower), interval0.well_defined);
}
IntervalSet arctanOpInterval(const Interval & interval) {
    return Interval(atan(interval.lower), atan(interval.upper), interval.well_defined);
}
IntervalSet absOpInterval(const Interval & interval) {
    if (interval.lower >= 0.0) {
        return interval;
    } else if (interval.upper <= 0.0) {
        return negOpInterval(interval);
    } else {
        double vmax = std::max(-interval.lower, interval.upper);
        return Interval(0.0, vmax, interval.well_defined);
    }
}
IntervalSet rootOpInterval(const Interval & interval, int n_root) {
    if (n_root == 0) {  // nan
        return IntervalSet();
    } else if (n_root < 0) {
        return IntervalSet(1.0) / rootOpInterval(interval, -n_root);
    } else if (n_root % 2 == 0) {  // even root
        if (interval.upper < 0.0) {
            return IntervalSet();
        } else if (interval.lower < 0.0) {
            return Interval(0.0, pow(interval.upper, 1.0 / n_root), false);
        } else {
            return Interval(pow(interval.lower, 1.0 / n_root), pow(interval.upper, 1.0 / n_root), interval.well_defined);
        }
    } else {
        if (interval.upper < 0.0) {  // odd root
            return Interval(-pow(-interval.lower, 1.0 / n_root), -pow(-interval.upper, 1.0 / n_root), interval.well_defined);
        } else if (interval.lower < 0.0) {
            return Interval(-pow(-interval.lower, 1.0 / n_root), pow(interval.upper, 1.0 / n_root), interval.well_defined);
        } else {
            return Interval(pow(interval.lower, 1.0 / n_root), pow(interval.upper, 1.0 / n_root), interval.well_defined);
        }
    }
}
IntervalSet sqrtOpInterval(const Interval & interval) {
    return rootOpInterval(interval, 2);
}
IntervalSet powOpInterval(const Interval & interval1, const Interval & interval2) {
    Interval interval0;
    if (interval1.upper < 0.0) {
        return IntervalSet();
    } else if (interval1.lower < 0.0) {
        interval0 = Interval(0.0, interval1.upper, false);
    } else {
        interval0 = interval1;
    }
    bool defined;
    if (interval0.lower == 0.0 && interval2.lower <= 0.0) {
        defined = false;
    } else {
        defined = interval0.well_defined && interval2.well_defined;
    }
    return Interval(
        std::min({pow(interval0.lower, interval2.lower), pow(interval0.upper, interval2.upper), pow(interval0.lower, interval2.upper), pow(interval0.upper, interval2.lower)}),
        std::max({pow(interval0.lower, interval2.lower), pow(interval0.upper, interval2.upper), pow(interval0.lower, interval2.upper), pow(interval0.upper, interval2.lower)}),
        defined
    );
}
IntervalSet powIntOpInterval(const Interval & interval, int exponent) {
    if (exponent == 0) {
        if (interval.upper < 0.0) {
            return IntervalSet();
        } else {
            return Interval(1.0, interval.well_defined && interval.lower > 0.0);
        }
    } else if (exponent < 0) {
        return IntervalSet(1.0) / powIntOpInterval(interval, -exponent);
    } else if (exponent == 1) {
        return interval;
    } else {  // TODO: currently calculate recursively, but may be slow when exponent is large
        int exponent_half = exponent / 2;
        return powIntOpInterval(interval, exponent_half) * powIntOpInterval(interval, exponent - exponent_half);
    }
}
IntervalSet dblPowOpInterval(double base, const Interval & interval) {
    if (base < 0.0) {
        return IntervalSet();
    } else if (base == 0.0) {
        if (interval.upper <= 0.0) {
            return IntervalSet();
        } else if (interval.lower <= 0.0) {
            return Interval(0.0, 0.0, false);
        } else {
            return Interval(0.0, 0.0, interval.well_defined);
        }
    } else {
        return Interval(exp(interval.lower * log(base)), exp(interval.upper * log(base)), interval.well_defined);
    }
}

// operators like double + IntervalSet
IntervalSet operator + (double a, const IntervalSet & b) {
    return IntervalSet(a) + b;
}
IntervalSet operator - (double a, const IntervalSet & b) {
    return IntervalSet(a) - b;
}
IntervalSet operator * (double a, const IntervalSet & b) {
    return IntervalSet(a) * b;
}
IntervalSet operator / (double a, const IntervalSet & b) {
    return IntervalSet(a) / b;
}

IntervalSet expOp(const IntervalSet & a) {
    return a.applyUnaryOp(expOpInterval);
}
IntervalSet lnOp(const IntervalSet & a) {
    return a.applyUnaryOp(lnOpInterval);
}
IntervalSet sinOp(const IntervalSet & a) {
    return a.applyUnaryOp(sinOpInterval);
}
IntervalSet cosOp(const IntervalSet & a) {
    return a.applyUnaryOp(cosOpInterval);
}
IntervalSet tanOp(const IntervalSet & a) {
    return sinOp(a) / cosOp(a);
}
IntervalSet cotOp(const IntervalSet & a) {
    return cosOp(a) / sinOp(a);
}
IntervalSet secOp(const IntervalSet & a) {
    return 1.0 / cosOp(a);
}
IntervalSet cscOp(const IntervalSet & a) {
    return 1.0 / sinOp(a);
}
IntervalSet arcsinOp(const IntervalSet & a) {
    return a.applyUnaryOp(arcsinOpInterval);
}
IntervalSet arccosOp(const IntervalSet & a) {
    return a.applyUnaryOp(arccosOpInterval);
}
IntervalSet arctanOp(const IntervalSet & a) {
    return a.applyUnaryOp(arctanOpInterval);
}
IntervalSet absOp(const IntervalSet & a) {
    return a.applyUnaryOp(absOpInterval);
}
IntervalSet rootOp(const IntervalSet & a, int n_root) {
    return a.applyUnaryOp([n_root](const Interval & b) { return rootOpInterval(b, n_root); });
}
IntervalSet sqrtOp(const IntervalSet & a) {
    return a.applyUnaryOp(sqrtOpInterval);
}
IntervalSet powOp(const IntervalSet & a, const IntervalSet & b) {
    return a.applyBinaryOp(b, powOpInterval);
}
IntervalSet powOp(const IntervalSet & a, double exponent) {
    return a.applyUnaryOp([exponent](const Interval & base) { return powOpInterval(base, Interval(exponent)); });
}
IntervalSet powOp(const IntervalSet & a, int exponent) {
    return a.applyUnaryOp([exponent](const Interval & base) { return powIntOpInterval(base, exponent); });
}
IntervalSet powOp(double base, const IntervalSet & interval) {
    return interval.applyUnaryOp([base](const Interval & exponent) { return dblPowOpInterval(base, exponent); });
}

CmpResult cmp_not(CmpResult a) {
    if (a == TRUE) {
        return FALSE;
    } else if (a == FALSE) {
        return TRUE;
    } else {
        return MAYBE;
    }
}

CmpResult CmpEqual::operator () (const Interval & interval1, const Interval & interval2) {
    if (interval1.lower > interval2.upper || interval1.upper < interval2.lower) {
        return FALSE;
    } else {
        return MAYBE;
    }
}
CmpResult CmpEqual::operator () (const IntervalSet & intv_set1, const IntervalSet & intv_set2) {
    if (intv_set1.isEmpty() || intv_set2.isEmpty()) {
        return FALSE;
    }
    if (intv_set1.getLowerBound() > intv_set2.getUpperBound() || intv_set1.getUpperBound() < intv_set2.getLowerBound()) {
        return FALSE;
    } else {
        return MAYBE;
    }
}
bool CmpEqual::verifySolution(const Interval & interval1, const Interval & interval2, Func2D func) {
    bool has_pos = false;
    bool has_neg = false;
    IntervalSet v[4];
    v[0] = func(interval1.lower, interval2.lower);
    v[1] = func(interval1.lower, interval2.upper);
    v[2] = func(interval1.upper, interval2.lower);
    v[3] = func(interval1.upper, interval2.upper);
    for (int i = 0; i < 4; i++) {
        if (v[i].isEmpty()) {
            continue;
        }
        double value = v[i].intervals[0].lower;
        if (value > 0.0) {
            has_pos = true;
        } else if (value < 0.0) {
            has_neg = true;
        } else {
            return true;
        }
    }
    return has_pos && has_neg;
}

CmpResult CmpLesser::operator () (const Interval & interval1, const Interval & interval2) {
    if (interval1.lower > interval2.upper) {
        return FALSE;
    } else if (interval1.upper < interval2.lower) {
        return TRUE;
    } else {
        return MAYBE;
    }
}
CmpResult CmpLesser::operator () (const IntervalSet & a, const IntervalSet & b) {
    if (a.isEmpty() || b.isEmpty()) {
        return FALSE;
    }
    if (a.getLowerBound() > b.getUpperBound()) {
        return FALSE;
    } else if (a.getUpperBound() < b.getLowerBound()) {
        return TRUE;
    } else {
        return MAYBE;
    }
}
bool CmpLesser::verifySolution(const Interval & interval1, const Interval & interval2, Func2D func) {
    IntervalSet v[4];
    v[0] = func(interval1.lower, interval2.lower);
    v[1] = func(interval1.lower, interval2.upper);
    v[2] = func(interval1.upper, interval2.lower);
    v[3] = func(interval1.upper, interval2.upper);
    for (int i = 0; i < 4; i++) {
        if (v[i].isEmpty()) {
            continue;
        }
        double value = v[i].intervals[0].lower;
        if (value < 0.0) {
            return true;
        }
    }
    return false;
}

CmpResult CmpGreater::operator () (const Interval & interval1, const Interval & interval2) {
    return CmpLesser()(interval2, interval1);
}
CmpResult CmpGreater::operator () (const IntervalSet & a, const IntervalSet & b) {
    return CmpLesser()(b, a);
}
bool CmpGreater::verifySolution(const Interval & interval1, const Interval & interval2, Func2D func) {
    IntervalSet v[4];
    v[0] = func(interval1.lower, interval2.lower);
    v[1] = func(interval1.lower, interval2.upper);
    v[2] = func(interval1.upper, interval2.lower);
    v[3] = func(interval1.upper, interval2.upper);
    for (int i = 0; i < 4; i++) {
        if (v[i].isEmpty()) {
            continue;
        }
        double value = v[i].intervals[0].lower;
        if (value > 0.0) {
            return true;
        }
    }
    return false;
}

CmpResult CmpLesserEq::operator () (const Interval & interval1, const Interval & interval2) {
    return cmp_not(CmpGreater()(interval1, interval2));
}
CmpResult CmpLesserEq::operator () (const IntervalSet & a, const IntervalSet & b) {
    return cmp_not(CmpGreater()(a, b));
}
bool CmpLesserEq::verifySolution(const Interval & interval1, const Interval & interval2, Func2D func) {
    IntervalSet v[4];
    v[0] = func(interval1.lower, interval2.lower);
    v[1] = func(interval1.lower, interval2.upper);
    v[2] = func(interval1.upper, interval2.lower);
    v[3] = func(interval1.upper, interval2.upper);
    for (int i = 0; i < 4; i++) {
        if (v[i].isEmpty()) {
            continue;
        }
        double value = v[i].intervals[0].lower;
        if (value <= 0.0) {
            return true;
        }
    }
    return false;
}

CmpResult CmpGreaterEq::operator () (const Interval & interval1, const Interval & interval2) {
    return cmp_not(CmpLesser()(interval1, interval2));
}
CmpResult CmpGreaterEq::operator () (const IntervalSet & a, const IntervalSet & b) {
    return cmp_not(CmpLesser()(a, b));
}
bool CmpGreaterEq::verifySolution(const Interval & interval1, const Interval & interval2, Func2D func) {
    IntervalSet v[4];
    v[0] = func(interval1.lower, interval2.lower);
    v[1] = func(interval1.lower, interval2.upper);
    v[2] = func(interval1.upper, interval2.lower);
    v[3] = func(interval1.upper, interval2.upper);
    for (int i = 0; i < 4; i++) {
        if (v[i].isEmpty()) {
            continue;
        }
        double value = v[i].intervals[0].lower;
        if (value >= 0.0) {
            return true;
        }
    }
    return false;
}

void Func2DPlotter::writeFile(const std::string & filepath) const {
    std::ofstream ofs(filepath, std::ios::binary);
    if (!ofs) {
        std::cerr << "cannot open file " << filepath << "!" << std::endl;
        return;
    }
    ofs.write(reinterpret_cast<const char *>(this->data), sizeof(int) * this->Nx * this->Ny);
    ofs.close();
}

void Func2DPlotter::plotRecur(int x1, int x2, int y1, int y2) {
    IntervalSet equation = this->func(this->sampleCoordX(x1, x2), this->sampleCoordY(y1, y2));
    CmpResult eq_cmp = (*(this->relation))(equation, 0.0);
    if (x2 - x1 <= 1 && y2 - y1 <= 1) {
        if (!equation.isWellDefined()) {
            return;
        } else if (eq_cmp == TRUE) {
            this->setPixel(x1, y1, 1);
        } else if (eq_cmp == MAYBE && this->verifySolution(x1, y1)) {
            this->setPixel(x1, y1, 1);
        }
    } else if (eq_cmp == TRUE && equation.isWellDefined()) {
        for (int x = x1; x < x2; x++) {
            for (int y = y1; y < y2; y++) {
                this->setPixel(x, y, 1);
            }
        }
    } else if (eq_cmp != FALSE) {
        int xm = (x1 + x2) / 2;
        int ym = (y1 + y2) / 2;
        this->plotRecur(x1, xm, y1, ym);
        this->plotRecur(x1, xm, ym, y2);
        this->plotRecur(xm, x2, y1, ym);
        this->plotRecur(xm, x2, ym, y2);
    }
}
