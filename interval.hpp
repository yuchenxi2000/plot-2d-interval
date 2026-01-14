//
// 2D Graphing by Interval Arithmetic, Implemented in C++
// References: https://zhuanlan.zhihu.com/p/1988701719459284501 and https://www.dgp.toronto.edu/public_user/mooncake/papers/SIGGRAPH2001_Tupper.pdf
// Copyright (c) 2026 YCX. Licensed under MIT License.
//

#ifndef interval_hpp
#define interval_hpp

#include <vector>
#include <algorithm>
#include <limits>
#include <functional>
#include <string>


enum CmpResult {
    TRUE, FALSE, MAYBE
};

class Interval {
public:
    double lower, upper;
    bool well_defined;
    Interval(double a, double b, bool defined) {
        lower = a;
        upper = b;
        well_defined = defined;
    }
    Interval(double a, double b) {
        lower = a;
        upper = b;
        well_defined = true;
    }
    Interval(double a) {
        lower = a;
        upper = a;
        well_defined = true;
    }
    Interval() {
        well_defined = true;
    }
    static Interval whole(bool defined = true) {
        return Interval(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), defined);
    }
    static Interval toPosInf(double a, bool defined = true) {
        return Interval(a, std::numeric_limits<double>::infinity(), defined);
    }
    static Interval toNegInf(double a, bool defined = true) {
        return Interval(-std::numeric_limits<double>::infinity(), a, defined);
    }
};

class IntervalSet;
typedef std::function<IntervalSet(const Interval &)> UnaryOp;
typedef std::function<IntervalSet(const Interval &, const Interval &)> BinaryOp;

// addition, subtraction, multiplication, division
IntervalSet addOpInterval(const Interval & a, const Interval & b);
IntervalSet subOpInterval(const Interval & a, const Interval & b);
IntervalSet negOpInterval(const Interval & a);
IntervalSet mulOpInterval(const Interval & a, const Interval & b);
IntervalSet divOpInterval(const Interval & a, const Interval & b);

class IntervalSet {
public:
    std::vector<Interval> intervals;
    IntervalSet() {}
    IntervalSet(const std::vector<Interval> & intvs) : intervals(intvs) {}
    IntervalSet(const Interval & interval) {
        intervals = std::vector<Interval>({interval});
    }
    IntervalSet(double a) {
        intervals = std::vector<Interval>({Interval(a)});
    }
    const std::vector<Interval> & getIntervals() const { return intervals; }

    bool isWellDefined() const;
    bool isEmpty() const {
        return intervals.empty();
    }
    static IntervalSet empty() {
        return IntervalSet(std::vector<Interval>());
    }
    static IntervalSet whole(bool defined = true) {
        return IntervalSet({Interval::whole(defined)});
    }

    IntervalSet applyUnaryOp(UnaryOp func) const;
    IntervalSet applyBinaryOp(const IntervalSet & other, BinaryOp func) const;

    double getLowerBound() const;
    double getUpperBound() const;
    
    // merge intervals that overlap
    IntervalSet mergeOverlapping() const;

    IntervalSet operator + (const IntervalSet & other) const {
        return this->applyBinaryOp(other, addOpInterval);
    }
    IntervalSet operator - (const IntervalSet & other) const {
        return this->applyBinaryOp(other, subOpInterval);
    }
    IntervalSet operator - () const {
        return this->applyUnaryOp(negOpInterval);
    }
    IntervalSet operator * (const IntervalSet & other) const {
        return this->applyBinaryOp(other, mulOpInterval);
    }
    IntervalSet operator / (const IntervalSet & other) const {
        return this->applyBinaryOp(other, divOpInterval);
    }
};

// operators like double + IntervalSet
IntervalSet operator + (double a, const IntervalSet & b);
IntervalSet operator - (double a, const IntervalSet & b);
IntervalSet operator * (double a, const IntervalSet & b);
IntervalSet operator / (double a, const IntervalSet & b);

// exponent, logrithm, (arc) triangular functions, absolute value, n-th root, power
IntervalSet expOp(const IntervalSet & a);
IntervalSet lnOp(const IntervalSet & a);
IntervalSet sinOp(const IntervalSet & a);
IntervalSet cosOp(const IntervalSet & a);
IntervalSet tanOp(const IntervalSet & a);
IntervalSet cotOp(const IntervalSet & a);
IntervalSet secOp(const IntervalSet & a);
IntervalSet cscOp(const IntervalSet & a);
IntervalSet arcsinOp(const IntervalSet & a);
IntervalSet arccosOp(const IntervalSet & a);
IntervalSet arctanOp(const IntervalSet & a);
IntervalSet absOp(const IntervalSet & a);
IntervalSet rootOp(const IntervalSet & a, int n_root);
IntervalSet sqrtOp(const IntervalSet & a);
IntervalSet powOp(const IntervalSet & a, const IntervalSet & b);
IntervalSet powOp(const IntervalSet & a, double exponent);
IntervalSet powOp(const IntervalSet & a, int exponent);
IntervalSet powOp(double base, const IntervalSet & interval);

CmpResult cmp_not(CmpResult a);
typedef std::function<IntervalSet(const IntervalSet &, const IntervalSet &)> Func2D;

class Cmp {
public:
    Cmp() {}
    virtual ~Cmp() {}
    virtual CmpResult operator () (const Interval & interval1, const Interval & interval2) = 0;
    virtual CmpResult operator () (const IntervalSet & interval1, const IntervalSet & interval2) = 0;
    virtual bool verifySolution(const Interval & interval1, const Interval & interval2, Func2D func) = 0;
};

class CmpEqual : public Cmp {
public:
    virtual CmpResult operator () (const Interval & interval1, const Interval & interval2);
    virtual CmpResult operator () (const IntervalSet & intv_set1, const IntervalSet & intv_set2);
    virtual bool verifySolution(const Interval & interval1, const Interval & interval2, Func2D func);
};

class CmpLesser : public Cmp {
public:
    virtual CmpResult operator () (const Interval & interval1, const Interval & interval2);
    virtual CmpResult operator () (const IntervalSet & a, const IntervalSet & b);
    virtual bool verifySolution(const Interval & interval1, const Interval & interval2, Func2D func);
};

class CmpGreater : public Cmp {
public:
    virtual CmpResult operator () (const Interval & interval1, const Interval & interval2);
    virtual CmpResult operator () (const IntervalSet & a, const IntervalSet & b);
    virtual bool verifySolution(const Interval & interval1, const Interval & interval2, Func2D func);
};

class CmpLesserEq : public Cmp {
public:
    virtual CmpResult operator () (const Interval & interval1, const Interval & interval2);
    virtual CmpResult operator () (const IntervalSet & a, const IntervalSet & b);
    virtual bool verifySolution(const Interval & interval1, const Interval & interval2, Func2D func);
};

class CmpGreaterEq : public Cmp {
public:
    virtual CmpResult operator () (const Interval & interval1, const Interval & interval2);
    virtual CmpResult operator () (const IntervalSet & a, const IntervalSet & b);
    virtual bool verifySolution(const Interval & interval1, const Interval & interval2, Func2D func);
};

class Func2DPlotter {
public:
    Func2D func;
    Cmp * relation;
    Interval interval_x;
    Interval interval_y;
    int Nx;
    int Ny;
    int * data;
    Func2DPlotter(Func2D func, Cmp * relation, const Interval & interval_x, const Interval & interval_y, int Nx, int Ny) {
        this->func = func;
        this->relation = relation;
        this->interval_x = interval_x;
        this->interval_y = interval_y;
        this->Nx = Nx;
        this->Ny = Ny;
        this->data = new int[Nx * Ny]();
    }

    ~Func2DPlotter() {
        delete [] this->data;
    }

    Interval sampleCoordX(int a, int b) const {
        double w1 = (a - 0.5) / (Nx - 1);
        double w2 = (b - 0.5) / (Nx - 1);
        return Interval(interval_x.lower + (interval_x.upper - interval_x.lower) * w1, interval_x.lower + (interval_x.upper - interval_x.lower) * w2);
    }

    Interval sampleCoordY(int a, int b) const {
        double w1 = (a - 0.5) / (Ny - 1);
        double w2 = (b - 0.5) / (Ny - 1);
        return Interval(interval_y.lower + (interval_y.upper - interval_y.lower) * w1, interval_y.lower + (interval_y.upper - interval_y.lower) * w2);
    }

    void setPixel(int x, int y, int value) {
        this->data[y + x * this->Ny] = value;
    }

    void writeFile(const std::string & filepath) const;

    bool verifySolution(int x, int y) const {
        Interval interval1 = this->sampleCoordX(x, x + 1);
        Interval interval2 = this->sampleCoordY(y, y + 1);
        return this->relation->verifySolution(interval1, interval2, this->func);
    }

    void plotRecur(int x1, int x2, int y1, int y2);

    void plot() {
        this->plotRecur(0, this->Nx, 0, this->Ny);
    }
};

#endif /* interval_hpp */
