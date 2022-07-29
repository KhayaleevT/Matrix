//
// Created by maent on 26.03.2020.
//
// Created by maent on 25.03.2020.
//

#ifndef TEMPLATES_POLYNOMIAL_MULTIPLY_H
#define TEMPLATES_POLYNOMIAL_MULTIPLY_H

#include <iostream>
#include <utility>
#include <cassert>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

using std::pair;
using std::vector;
using std::complex;
typedef long long LL;
typedef complex<long double> complex_d;

void fft(vector<complex_d> &a, bool invert, complex_d root) {
    int n = a.size();
    if (n == 1) {
        return;
    }
    vector<complex_d> a0(n / 2);
    vector<complex_d> a1(n / 2);
    for (int i = 0, j = 0; i < n; i += 2, j++) {
        a0[j] = a[i];
        a1[j] = a[i + 1];
    }
    fft(a0, invert, root * root);
    fft(a1, invert, root * root);
    complex_d r = 1;
    complex_d complex_2 = 2;
    for (int i = 0; i < n / 2; i++) {
        a[i] = a0[i] + r * a1[i];
        a[i + n / 2] = a0[i] - r * a1[i];
        if (invert) {
            a[i] /= complex_2;
            a[i + n / 2] /= complex_2;
        }
        r *= root;
    }
}

pair<int, int> two_power_ceil(int n) {
    int a = 0, pow = 1;
    while (pow < n) {
        pow *= 2;
        a++;
    }
    return {a, pow};
}


vector<long long> Multiply_Polynoms(const vector<long long> &a, const vector<long long> &b) {
    int n = int(a.size()) + int(b.size()) - 1;
    pair<int, int> deg = two_power_ceil(n);
    vector<complex_d> ca(deg.second);
    vector<complex_d> cb(deg.second);
    for (size_t i = 0; i < a.size(); i++) {
        ca[i] = complex_d(a[i], 0);
    }
    for (size_t i = 0; i < b.size(); i++) {
        cb[i] = complex_d(b[i], 0);
    }
    double ang(2 * M_PI / deg.second);
    complex_d primitive_root(cos(ang), sin(ang));
    fft(ca, false, primitive_root);
    fft(cb, false, primitive_root);
    for (size_t i = 0; i < ca.size(); i++) {
        ca[i] *= cb[i];
    }
    complex_d inv_r = complex_d(1) / primitive_root;
    fft(ca, true, inv_r);
    vector<long long> res(ca.size());
    for (size_t i = 0; i < res.size(); i++) {
        res[i] = static_cast<long long>(ca[i].real() + 0.5);
    }
    while (res.back() == 0) {
        res.pop_back();
    }
    return res;
}

const long long N_BASE = int(1e6);

void multiply(vector<long long> &a, vector<long long> &b, vector<long long> &ans) {
    vector<long long> C = Multiply_Polynoms(a, b);
    C.resize(a.size() + b.size(), 0);
    for (unsigned i = 0; i + 1 < C.size(); i++) {
        int digit = C[i] % N_BASE;
        C[i + 1] += ((C[i] - LL(digit)) / N_BASE);
        C[i] = digit;
    }
    for (int i = int(C.size()) - 1; (C[i] == 0) && (i >= 0); i--) {
        C.pop_back();
    }
    ans = C;
}


#endif //TEMPLATES_POLYNOMIAL_MULTIPLY_H

#pragma once

#include <iostream>
#include<deque>
#include<cassert>

using std::string;
using std::istream;
using std::ostream;

#ifndef RATIONAL
#define RATIONAL

const int NSYS_LOG = 6;
#ifndef MY_COMPARE
#define MY_COMPARE
enum myComparator {
    GREATER = 1, LESS = -1, EQUAL = 0
};
#endif
const int NEGATIVE = 0;
const int NON_NEGATIVE = 1;

class BigInteger {
private:
    bool isSmall() const {
        return _num.size() <= 3;
    }

public:
    void change_sign() {
        _sign = !_sign;
    }

    bool _IsNull() const {
        return (_num.size() == 1) && (_num[0] == 0);
    }

    bool isNegative() const {
        return !_sign;
    }

    bool IsUnit() const {
        return (_num.size() == 1) && (_num[0] == 1);
    }

public:
    friend istream &operator>>(istream &in, BigInteger &x) {
        string s;
        in >> s;
        x = BigInteger(s);
        return in;
    }

    friend ostream &operator<<(ostream &out, const BigInteger &x) {
        if (x._sign == NEGATIVE)out << "-";
        for (int i = int(x._num.size()) - 1; i >= 0; i--) {
            int deg = 10;
            for (int j = 1; j < NSYS_LOG; j++) {
                if (x._num[i] < deg && (i != int(x._num.size()) - 1)) {
                    out << "0";
                }
                deg *= 10;
            }
            out << x._num[i];
        }
        return out;
    }

    BigInteger(const string &s) {
        _sign = NON_NEGATIVE;
        if (s[0] == '-') {
            _sign = NEGATIVE;
        }
        for (int i = int(s.length()) - 1; i >= int(!_sign); i -= NSYS_LOG) {
            int resid = 0;
            int deg = 1;
            for (int j = 0; (j < NSYS_LOG) && (i - j >= int(!_sign)); j++) {
                resid += (deg * (s[i - j] - '0'));
                deg *= 10;
            }
            _num.push_back(resid);
        }
    }


    BigInteger() {
        _sign = NON_NEGATIVE;
        _num = {0};
    }

    BigInteger(int x) {
        if (x == 0) {
            _num.push_back(0);
            _sign = NON_NEGATIVE;
            return;
        }
        _sign = NON_NEGATIVE;
        if (x < 0) {
            _sign = NEGATIVE;
            x = -x;
        }
        while (x > 0) {
            _num.push_back(x % N_BASE);
            x /= N_BASE;
        }
    }

    size_t amount_of_digits() {
        return _num.size();
    }

    long long toLL(int l, int r) {
        long long deg = 1;
        long long ans = 0;
        for (int i = l; i <= r; i++) {
            ans += (_num[i] * deg);
            deg *= N_BASE;
        }
        return ans;
    }

    explicit operator bool() const {
        return !_IsNull();
    }

    string toString(int l = 0, int r = -10, bool WithSign = true, bool IfFrac = false) const {
        if (r == -10) {
            r = int(_num.size()) - 1;
        }
        if (l > r) {
            if (WithSign && !(*this)._sign)return "-0";
            return "0";
        }
        string a;
        if (WithSign && !(*this)._sign)a.push_back('-');
        int i = r;
        while (i >= int(_num.size())) {
            for (int j = 0; j < NSYS_LOG; j++) {
                a.push_back('0');
            }
            i--;
        }
        for (; i >= l; i--) {
            int deg = 10;
            for (int j = 1; j < NSYS_LOG; j++) {
                if ((*this)._num[i] < deg && (i != int((*this)._num.size()) - 1 || IfFrac)) {
                    a.push_back('0');
                }
                deg *= 10;
            }
            a += std::to_string(_num[i]);

        }
        return a;
    }


    friend bool operator>(const BigInteger &a, const BigInteger &b) {
        return (_cmp(a, b) == GREATER);
    }

    friend bool operator>=(const BigInteger &a, const BigInteger &b) {
        return (_cmp(a, b) != LESS);
    }

    friend bool operator<(const BigInteger &a, const BigInteger &b) {
        return (_cmp(a, b) == LESS);
    }

    friend bool operator<=(const BigInteger &a, const BigInteger &b) {
        return (_cmp(a, b) != GREATER);
    }

    friend bool operator==(const BigInteger &a, const BigInteger &b) {
        return (_cmp(a, b) == EQUAL);
    }

    friend bool operator!=(const BigInteger &a, const BigInteger &b) {
        return (_cmp(a, b) != EQUAL);
    }

    BigInteger operator-() const {
        BigInteger x = *this;
        if (x == 0) {
            return x;
        }
        x._sign ^= 1;
        return x;
    }

    BigInteger &operator+=(const BigInteger &other) {
        _sum(*this, other, *this);
        return (*this);
    }

    BigInteger &operator-=(const BigInteger &other) {
        if (other == (BigInteger) 0) {
            return (*this);
        }
        _sum(*this, other, *this, true);
        return (*this);
    }

    friend BigInteger operator+(const BigInteger &a, const BigInteger &b) {
        BigInteger ans = a;
        ans += b;
        return ans;
    }

    friend BigInteger operator-(const BigInteger &a, const BigInteger &b) {
        BigInteger ans = a;
        ans -= b;
        return ans;
    }

    BigInteger &operator++() {
        return *this += 1;
    }

    BigInteger operator++(int) {
        ++(*this);
        return *this - 1;
    }

    BigInteger &operator--() {
        return *this -= 1;
    }

    BigInteger operator--(int) {
        --(*this);
        return *this + 1;
    }

    friend void multiply_by_fft(const BigInteger &a, const BigInteger &b, BigInteger &ans) {
        vector<long long> A(a._num.size());
        vector<long long> B(b._num.size());
        for (unsigned i = 0; i < A.size(); i++) {
            A[i] = a._num[i];
        }
        for (unsigned i = 0; i < B.size(); i++) {
            B[i] = b._num[i];
        }
        vector<long long> ANS;
        multiply(A, B, ANS);
        ans._num = {};
        for (unsigned i = 0; i < ANS.size(); i++) {
            ans._num.push_back(ANS[i]);
        }
    }

    friend BigInteger operator*(const BigInteger &a, const BigInteger &b) {
        BigInteger ans = 0;
        if (a._IsNull() || b._IsNull())return ans;
        if (a.isSmall()) {
            _positive_multiplication(b, a, ans);
        } else if (b.isSmall()) {
            _positive_multiplication(a, b, ans);
        } else {
            multiply_by_fft(a, b, ans);
        }
        if (a._sign != b._sign && (a != 0) && (b != 0)) {
            ans._sign = NEGATIVE;
            return ans;
        }
        return ans;
    }

    BigInteger &operator*=(const BigInteger &b) {
        (*this) = (*this) * b;
        return *this;
    }

    BigInteger &operator/=(const BigInteger &b) {
        if (_IsNull() || b == 1) {
            return *this;
        }
        if (_sign == b._sign) {
            (*this)._sign = NON_NEGATIVE;
            *this = _div_for_absolutes(*this, b).first;
            return *this;
        }
        *this = _div_for_absolutes(*this, b).first;
        if (!_IsNull()) {
            (*this)._sign = NEGATIVE;
        }
        return *this;
    }

    friend BigInteger operator/(const BigInteger &a, const BigInteger &b) {
        BigInteger ans = a;
        ans /= b;
        return ans;
    }

    BigInteger &operator%=(const BigInteger &b) {
        int Original_Sign = _sign;
        (*this) = _div_for_absolutes(*this, b).second;
        if (!_IsNull() && Original_Sign != b._sign) {
            return *this *= -1;
        }
        return *this;
    }

    friend BigInteger operator%(const BigInteger &a, const BigInteger &b) {
        BigInteger ans = a;
        return ans %= b;
    }

    BigInteger &addZeros(const int &amount) {
        if (*this == 0)return *this;
        for (int i = 0; i < amount; i++) {
            _num.push_front(0);
        }
        return *this;
    }

    friend void swap(BigInteger &a, BigInteger &b) {
        std::swap(a._num, b._num);
        std::swap(a._sign, b._sign);
    }

    friend BigInteger abs(const BigInteger &a) {
        BigInteger absolute = a;
        absolute._sign = NON_NEGATIVE;
        return absolute;
    }

    friend BigInteger gcd(const BigInteger &a, const BigInteger &b) {
        BigInteger a1 = abs(a);
        BigInteger b1 = abs(b);
        int _comp = _cmp(a1, b1);
        if (_comp == LESS) {
            swap(a1, b1);
        }
        if (_comp == EQUAL) {
            return a1;
        }
        if (b1.IsUnit()) {
            return 1;
        }
        while (!a1._IsNull() && !b1._IsNull()) {
            a1 = _div_for_absolutes(a1, b1).second;
            if (a1 * BigInteger(2) > b1) {
                a1 = b1 - a1;
            }
            swap(a1, b1);
        }
        return a1;
    }

private:
    std::deque<int> _num;
    bool _sign;

    //sign='>=0'


    static int _cmp_for_absolutes_with_zeros_added_to_second(const BigInteger &a, const BigInteger &b,
                                                             size_t zeros_added_to_b = 0) {
        if (b._IsNull() && a._IsNull())return EQUAL;
        if (b._IsNull())return GREATER;
        int ans = EQUAL;
        if (a._num.size() > b._num.size() + zeros_added_to_b) {
            ans = GREATER;
            return ans;
        }
        if (a._num.size() < b._num.size() + zeros_added_to_b) {
            ans = LESS;
            return ans;
        }
        for (int i = 1; i <= int(b._num.size()); i++) {
            if (a._num[a._num.size() - i] > b._num[b._num.size() - i]) {
                ans = GREATER;
                return ans;
            }
            if (b._num[b._num.size() - i] > a._num[a._num.size() - i]) {
                ans = LESS;
                return ans;
            }
        }
        for (int i = int(b._num.size()) + 1; i <= int(a._num.size()); i++) {
            if (a._num[a._num.size() - i] > 0) {
                ans = GREATER;
                return ans;
            }
        }
        return ans;
    }

    static int _cmp(const BigInteger &a, const BigInteger &b, bool ChangeFirstSign = false,
                    bool ChangeSecondSign = false) {
        if ((a._sign ^ ChangeFirstSign) > (b._sign ^ ChangeSecondSign)) {
            return GREATER;
        }
        if ((a._sign ^ ChangeFirstSign) < (b._sign ^ ChangeSecondSign)) {
            return LESS;
        }
        int ans = _cmp_for_absolutes_with_zeros_added_to_second(a, b);
        if (!(a._sign ^ ChangeFirstSign)) {
            ans = -ans;
        }
        return ans;
    }


    //(|a|+|b|) or |a|+(-|b|) only if |a|>=|b|
    static void _non_negative_Sum_or_Sub_for_absolutes(const BigInteger &a, const BigInteger &b, BigInteger &ans,
                                                       bool IfSub = false) {
        int resid = 0;
        int factor = 1;
        if (IfSub)factor = -1;
        size_t num_size = std::max(a._num.size(), b._num.size());
        for (size_t i = 0; (i < num_size) || (resid); i++) {
            int add = resid;
            if (i < a._num.size()) {
                add += a._num[i];
            }
            if (i < b._num.size()) {
                add += ((factor) * b._num[i]);
            }
            if (add >= N_BASE) {
                resid = 1;
                add -= N_BASE;
            } else if (add < 0) {
                resid = -1;
                add += N_BASE;
            } else {
                resid = 0;
            }
            if (ans._num.size() > i) {
                ans._num[i] = add;
            } else {
                ans._num.push_back(add);
            }
        }
        int s = ans._num.size();
        for (int i = s - 1; i >= 1; i--) {
            if (ans._num[i] != 0)break;
            ans._num.pop_back();
        }
    }

    //(a+b) or a+(-b)
    static void _sum(const BigInteger &a, const BigInteger &b, BigInteger &ans, bool ChangeSecondSign = false) {
        if (a._sign == (b._sign ^ ChangeSecondSign)) {
            _non_negative_Sum_or_Sub_for_absolutes(a, b, ans);
            ans._sign = a._sign;
        }
        if (a._sign && !(b._sign ^ ChangeSecondSign)) {
            if (_cmp(a, b, false, (true ^ ChangeSecondSign)) != -1) {
                ans._sign = NON_NEGATIVE;
                _non_negative_Sum_or_Sub_for_absolutes(a, b, ans, true);
            } else {
                _non_negative_Sum_or_Sub_for_absolutes(b, a, ans, true);
                ans._sign = NEGATIVE;
            }
        }
        if (!a._sign && (b._sign ^ ChangeSecondSign)) {
            if (_cmp(a, b, true, (false ^ ChangeSecondSign)) == 1) {
                _non_negative_Sum_or_Sub_for_absolutes(a, b, ans, true);
                ans._sign = NEGATIVE;
            } else {
                _non_negative_Sum_or_Sub_for_absolutes(b, a, ans, true);
                ans._sign = NON_NEGATIVE;
            }
        }
    }

    BigInteger abs_multiplied_by_digit(int x) const {
        if (x == 0)return 0;
        BigInteger ans = *this;
        ans._sign = NON_NEGATIVE;
        long long resid = 0;
        for (size_t i = 0; i < _num.size(); i++) {
            long long digit = _num[i];
            digit *= x;
            digit += resid;
            ans._num[i] = digit % N_BASE;
            resid = digit / N_BASE;
        }
        if (resid != 0) {
            ans._num.push_back(resid);
        }
        return ans;
    }

    static void _positive_multiplication(const BigInteger &a, const BigInteger &b, BigInteger &ans) {
        ans = 0;
        for (int i = int(b._num.size()) - 1; i >= 0; i--) {
            ans += ((a.abs_multiplied_by_digit(b._num[i])).addZeros(i));
        }
    }

    //ratio>1, a and b non_negative
    static int _NSYS_DEG_LOWER_BOUND_FOR_RATIO(const BigInteger &a, const BigInteger &b) {
        int NSYS_DEG_LOWER_BOUND = int(a._num.size()) - int(b._num.size());
        if (NSYS_DEG_LOWER_BOUND != 0) {
            if (_cmp_for_absolutes_with_zeros_added_to_second(a, b, NSYS_DEG_LOWER_BOUND) ==
                LESS)
                NSYS_DEG_LOWER_BOUND--;
        }
        return NSYS_DEG_LOWER_BOUND;
    }

    static int _FIND_FIRST_DIGIT_IN_DIV(const BigInteger &c, const BigInteger &b, int amount_of_digits_in_div) {
        int l = -1;
        int r = N_BASE;
        while (r - l > 1) {
            int m = (r + l) / 2;
            if (_cmp_for_absolutes_with_zeros_added_to_second(c, b.abs_multiplied_by_digit(m),
                                                              amount_of_digits_in_div - 1) != LESS) {
                l = m;
            } else {
                r = m;
            }
        }
        return l;
    }

    static std::pair<BigInteger, BigInteger> _div_for_absolutes(const BigInteger &a, const BigInteger &b) {
        BigInteger c = a;
        c._sign = NON_NEGATIVE;
        if (_cmp_for_absolutes_with_zeros_added_to_second(c, b, 0) == LESS) {
            return {0, c};
        }
        int NSYS_DEG_LOWER_BOUND = _NSYS_DEG_LOWER_BOUND_FOR_RATIO(c, b);
        std::pair<BigInteger, BigInteger> ans = {0, 0};
        ans.first._num.pop_back();
        for (int i = NSYS_DEG_LOWER_BOUND; i >= 0; i--) {
            int l = _FIND_FIRST_DIGIT_IN_DIV(c, b, i + 1);
            ans.first._num.push_front(l);
            c -= b.abs_multiplied_by_digit(l).addZeros(i);
        }
        ans.second = c;
        return ans;
    }
};


class Rational {
public:
    bool IsUnit() const {
        return _Numerator.IsUnit() && _Denominator.IsUnit();
    }

    bool IsZero() const {
        return _Numerator._IsNull();
    }

    friend istream &operator>>(istream &in, Rational &x) {
        string s;
        in >> s;
        x = BigInteger(s);
        return in;
    }

    Rational() {
        _Numerator = 0;
        _Denominator = 1;
    }

    Rational(const BigInteger &a) {
        _Numerator = a;
        _Denominator = 1;
    }

    Rational(int a) {
        _Numerator = a;
        _Denominator = 1;
    }

    friend bool operator==(const Rational &a, const Rational &b) {
        return (a._Numerator == b._Numerator) && (b._Denominator == a._Denominator);
    }

    friend bool operator!=(const Rational &a, const Rational &b) {
        return (a._Numerator != b._Numerator) || (b._Denominator != a._Denominator);
    }

    friend bool operator>=(const Rational &a, const Rational &b) {
        return (a._Numerator * b._Denominator >= a._Denominator * b._Numerator);
    }

    friend bool operator<=(const Rational &a, const Rational &b) {
        return (a._Numerator * b._Denominator <= a._Denominator * b._Numerator);
    }

    friend bool operator>(const Rational &a, const Rational &b) {
        return (a._Numerator * b._Denominator > a._Denominator * b._Numerator);
    }

    friend bool operator<(const Rational &a, const Rational &b) {
        return (a._Numerator * b._Denominator < a._Denominator * b._Numerator);
    }

    void reduce() {
        if (_Numerator.IsUnit() || _Denominator.IsUnit()) {
            return;
        }
        if (_Numerator._IsNull()) {
            _Denominator = 1;
            return;
        }
        BigInteger GCD = gcd(_Numerator, _Denominator);
        if (GCD.IsUnit()) {
            return;
        }
        _Numerator /= GCD;
        _Denominator /= GCD;
    }


    Rational &operator+=(const Rational &b) {
        Rational b1 = b;
        BigInteger GCD = gcd(_Denominator, b._Denominator);
        if (!GCD.IsUnit()) {
            _Denominator /= GCD;
            b1._Denominator /= GCD;
        }
        _Numerator = _Numerator * b1._Denominator + _Denominator * b1._Numerator;
        BigInteger save_denom = _Denominator;
        _Denominator = GCD;
        reduce();
        if (!_Numerator._IsNull()) {
            _Denominator *= (save_denom * b1._Denominator);
        }
        return *this;
    }

    friend Rational operator+(const Rational &a, const Rational &b) {
        Rational ans = a;
        ans += b;
        return ans;
    }

    Rational &operator-=(const Rational &b) {
        *this += (-b);
        return *this;
    }

    friend Rational operator-(const Rational &a, const Rational &b) {
        Rational ans = a;
        ans -= b;
        return ans;
    }

    Rational operator-() const {
        Rational ans = (*this);
        ans._Numerator.change_sign();
        return ans;
    }

    Rational &operator++() {
        return *this += 1;
    }

    Rational operator++(int) {
        ++(*this);
        return *this - Rational(1);
    }

    Rational &operator--() {
        return *this -= 1;
    }

    Rational operator--(int) {
        --(*this);
        return *this + Rational(1);
    }

    Rational &operator*=(const Rational &b) {
        Rational b1 = b;
        swap(b1._Denominator, _Denominator);
        reduce();
        b1.reduce();
        _Numerator *= b1._Numerator;
        _Denominator *= b1._Denominator;
        return *this;
    }

    friend Rational operator*(const Rational &a, const Rational &b) {
        Rational ans = a;
        ans *= b;
        return ans;
    }

    Rational inverted() const {
        assert((bool) _Numerator);
        Rational inv = *this;
        std::swap(inv._Numerator, inv._Denominator);
        if (inv._Denominator.isNegative()) {
            inv._Numerator.change_sign();
            inv._Denominator.change_sign();
        }
        return inv;
    }

    Rational &operator/=(const Rational &b) {
        *this *= b.inverted();
        return *this;
    }

    friend Rational operator/(const Rational &a, const Rational &b) {
        Rational ans = a;
        ans /= b;
        return ans;
    }

    string toString() {
        if (_Denominator == 1)return _Numerator.toString();
        string ans = _Numerator.toString();
        (ans += "/") += _Denominator.toString();
        return ans;
    }

    string asDecimal(size_t precision = 0) {
        BigInteger x = _Numerator.addZeros(precision);
        x /= _Denominator;
        string fractional_part = x.toString(0, int(precision) - 1, false, true);
        string floor = x.toString(int(precision), int(x.amount_of_digits()) - 1);
        return floor + "." + fractional_part;
    }


    explicit operator double() {
        return std::stod(asDecimal(50));
    }

private:
    BigInteger _Numerator, _Denominator;
};
//
// Created by maent on 15.03.2020.
//

#ifndef TEMPLATES_MATRIX_H
#define TEMPLATES_MATRIX_H

#include <iostream>
#include <vector>
#include<cassert>
#include <string>
#include <algorithm>

using std::pair;
using std::vector;
typedef long long LL;

template<unsigned n>
class Finite;


pair<long long, long long> linear_gcd_decomposition(long long a, long long b) {
    assert(a >= 0 && b >= 0);
    if (b == 0) {
        return {1, 0};
    }
    if (a == 0) {
        return {0, 1};
    }
    bool ifswap = false;
    if (a < b) {
        ifswap = true;
        std::swap(a, b);
    }
    long long k = a / b;
    pair<long long, long long> prev = linear_gcd_decomposition(a - k * b, b);
    pair<long long, long long> ans = {prev.first, prev.second - prev.first * k};
    if (ifswap) {
        std::swap(ans.first, ans.second);
    }
    return ans;
}

template<unsigned N>
class IsPrime {
private:
    template<unsigned N1, unsigned i>
    struct CheckDiv {
        static const bool check = CheckDiv<N, (N1 >= i * i) ? ((N1 % i == 0) ? 0 : i + 1) : 1>::check;
    };

    template<unsigned N1>
    struct CheckDiv<N1, 1> {
        static const bool check = true;
    };

    template<unsigned N1>
    struct CheckDiv<N1, 0> {
        static const bool check = false;
    };

public:
    static const bool value = CheckDiv<N, 2>::check;
};

template<>
class IsPrime<1> {
public:
    static const bool value = false;
};

template<unsigned n>
class upper_degree_of_two_bound {
private:
    template<unsigned N1, unsigned _pow>
    struct check_pow {
        static const unsigned pow = (N1 <= _pow) ? (_pow) : (check_pow<N1, _pow * 2>::pow);
    };
    template<unsigned N1>
    struct check_pow<N1, 2048> {
        static const unsigned pow = 2048;
    };
public:
    static const unsigned num = check_pow<n, 1>::pow;

};

template<bool u>
class My_Static_Assert {
public:
    static void assertion() {};
};

template<>
class My_Static_Assert<false> {
};

template<unsigned n>

class Finite {
protected:
    long long a = 0;
public:
    Finite() = default;

    Finite(long long a1) {
        a = ((a1 % n) + n) % n;
    }

    Finite &operator=(long long a1) {
        a = ((a1 % n) + n) % n;
        return *this;
    }

    Finite &operator*=(const Finite<n> &b) {
        a *= (b.a);
        a %= n;
        return *this;
    }

    Finite &operator+=(const Finite<n> &b) {
        a += b.a;
        a %= n;
        return *this;
    }

    Finite &operator++() {
        *this += Finite<n>(1);
        return *this;
    }

    Finite &operator-=(const Finite<n> &b) {
        a += (n - b.a);
        a %= n;
        return *this;
    }


    bool operator==(const Finite<n> &b) const {
        return a == b.a;
    }

    bool operator!=(const Finite<n> &b) const {
        return a != b.a;
    }

    explicit operator long long() const {
        return a;
    }

    Finite<n> inverse() const {
        My_Static_Assert<IsPrime<n>::value>::assertion();
        assert(a != 0);
        pair<long long, long long> inv = linear_gcd_decomposition(a, n);
        return Finite<n>(inv.first);
    }

    Finite<n> &operator/=(const Finite<n> &b) {
        My_Static_Assert<IsPrime<n>::value>::assertion();
        (*this) *= b.inverse();
        return *this;
    }
};

template<unsigned n>
const Finite<n> operator+(const Finite<n> &a, const Finite<n> &b) {
    Finite<n> c = a;
    c += b;
    return c;
}

template<unsigned n>
const Finite<n> operator*(const Finite<n> &a, const Finite<n> &b) {
    Finite<n> c = a;
    c *= b;
    return c;
}

template<unsigned n>
const Finite<n> operator-(const Finite<n> &a, const Finite<n> &b) {
    Finite<n> c = a;
    c -= b;
    return c;
}

template<unsigned n>
const Finite<n> operator/(const Finite<n> &a, const Finite<n> &b) {
    Finite<n> c = a;
    c /= b;
    return c;
}

template<typename T>
T Is_Field() {
    return T(1) / T(1);
}

template<unsigned M, unsigned N>
class static_max {
public:
    static const unsigned value = (M > N) ? M : N;
};

template<typename Field=Rational>
class MentallyHealthyMatrix;

template<unsigned M, unsigned N, typename Field=Rational>
class Matrix {
private:
    class Matrix_Row {
        friend Matrix<M, N, Field>;
    public:
        Field A[N];

        Matrix_Row() {
            for (unsigned i = 0; i < N; i++) {
                A[i] = Field(0);
            }
        };

        Matrix_Row(std::initializer_list<Field> a1) {
            int i = 0;
            for (auto c:a1) {
                A[i] = c;
                i++;
            }
        };

        Field &operator[](unsigned i) {
            return A[i];
        };

        const Field &operator[](unsigned i) const {
            return A[i];
        };
    };

    Matrix_Row _a[M];

    static pair<Field, int>
    _Gauss_operation(Matrix<M, N, Field> &A, unsigned starting_row, unsigned starting_column,
                     Matrix<M, N, Field> *Additional = nullptr) {
        int non_zero_row = -1;
        Field Mult = Field(0);
        Field CHANGE_MULTIPLY = Field(1);
        for (unsigned i = starting_row; i < M; i++) {
            if (non_zero_row == -1 && (A._a[i][starting_column] != Field(0))) {
                non_zero_row = (int) i;
                Mult = A._a[i][starting_column];
                if (i != starting_row) {
                    A.row_swap(i, starting_row);
                    if (Additional != nullptr) {
                        Additional->row_swap(i, starting_row);
                    }
                    CHANGE_MULTIPLY *= Field(-1);
                }
                continue;
            }
            if (non_zero_row != -1) {
                Field multiple = Field(-1) * (A._a[i][starting_column] / Mult);
                A.add_with_multiply(i, starting_row, multiple);
                A._a[i][starting_column] = Field(0);
                if (Additional != nullptr) {
                    Additional->add_with_multiply(i, starting_row, multiple);
                }
            }
        }
        return {CHANGE_MULTIPLY, non_zero_row};
    }
    /*





                ````````
            .///:-..-:/+o++++:.
         `-+/.``-:////oosyyhhhyo//:.`
       `-oy..:++:.``  ....````````-/+:.
      .osdo+:.`                     `-++.
    `/yyy:`                            -o/`
   `oNd:`                               `:o-
  -oom:                                   .+/`
  s/d+   `.--..`                .://///:.` `:o-                     `....`
 .hs/`  .yo:--://-`           -+s/----:/oo-.`.+/..                  yy+++.  ``            ``
.y/`   /yoo+//+++ss/`      `+s+/---------:+s.  -+s`                 do    :o+os.`y-  /s :soo+
y/    -o.   `/yss-`+/`     `s-   -ss/`    `/+   `h/`                do    .-::m/ oh`.m:.N:  .
/    `s.    `mNNm:  /+`    .o    oNMN/    ``y`   :h-                do    hs-:m+ `hoyo .N:  `
`    `s.    `:++:   `o/    `s`   `/+/.     `y`   `o`-.              s/    os++y/  .mh`  /yoo+
      -o.           .s-     :o.           .o-     o`:.                            -d.
       ./+:.`````.:++.       ./+-`     ``/+-      :/--                            ..
         `.://///:-`           `:////////-         o--
                                    `              /+-
      `.                                           :N`
     :oy`                                          oM.
    `-`y                                           oM.
       s-                                          sN.
       :s                                         /oy`
       `y.                                       .y:d:
        -o`                                      s/+o-:`
         ++                                     .y:+` ./`
         `:o.                     `             os:`   `/
           .++-`                 :hs-         `-o.      --
             .:+/:.`         `.-:+:``       `.::`       `:`
               ``.:///////////:-.``       `-:-`       `-omo`
                      ````             `-::.`      `-+dNMMN+
`                                  `.-::.`      `-ohNMMMMMMN.
o`                            `..-::-.`      .:odNMMMMMMMMMMo
/:` ```      ``      ```....::--..       `./sdNMMMMMMMMMMMMMd.
.ohhhdhhyysoshdyoyyhhddddds:`         .:ohmMMMMMMMMMMMMMMMMMM+
yMMMMMMMMMMMMMMMMMMMMMMMMMMNo`    `./ymNMMMMMMMMMMMMMMMMMMMMMd`
NMMMMMMMMMMMMMMMMMMMMMMMMMMMN. `./ymMMMMMMMMMMMMMMMMMMMMMMMMMM:
MMMMMMMMMNmmhdddo+oshmNMMMMMh`-smMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
MMMMMMdo/-.+ ``/-`  `.hmhddy:oNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
MMMMMMNho:+-`  /    `/NMh/-:dMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
MMMMMMMMMMNd/` /`  `oNMMMMNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
MMMMMMMMMMMNd/ /-``sNNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
MMMMMMMMMMMm` `/  `-+NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
MMMMMMMMMMMN/  /-``/NMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM+
hhhhhhhhhhhho  :``:yhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh:


     */
public:
    Matrix() {
        Is_Field<Field>();
        for (unsigned i = 0; i < M; i++) {
            for (unsigned j = 0; j < N; j++) {
                _a[i][j] = Field(0);
            }
        }
    }

    Matrix(Field x) : Matrix() {
        Is_Field<Field>();
        for (unsigned i = 0; i < std::min(M, N); i++) {
            _a[i][i] = x;
        }
    }

    Matrix(std::initializer_list<std::initializer_list<Field>> a1) {
        Is_Field<Field>();
        int i = 0;
        for (auto c:a1) {
            _a[i] = Matrix_Row(c);
            i++;
        }
    }

    Matrix_Row &operator[](unsigned i) {
        return _a[i];
    }

    const Matrix_Row &operator[](unsigned i) const {
        return _a[i];
    }

    Matrix &operator+=(const Matrix<M, N, Field> &B) {
        for (unsigned i = 0; i < M; i++) {
            for (unsigned j = 0; j < N; j++) {
                _a[i][j] += B[i][j];
            }
        }
        return *this;
    }

    Matrix &operator-=(const Matrix<M, N, Field> &B) {
        for (unsigned i = 0; i < M; i++) {
            for (unsigned j = 0; j < N; j++) {
                _a[i][j] -= B[i][j];
            }
        }
        return *this;
    }

    Matrix &operator*=(Field b) {
        for (unsigned i = 0; i < M; i++) {
            for (unsigned j = 0; j < N; j++) {
                _a[i][j] *= b;
            }
        }
        return *this;
    }


    Matrix<N, N, Field> &operator*=(const Matrix<N, N, Field> &b) {
        *this = classic_multiply(b);
        return *this;
    }

    void row_swap(unsigned i, unsigned j) {
        if (j == i)return;
        for (unsigned i1 = 0; i1 < N; i1++) {
            std::swap(_a[i][i1], _a[j][i1]);
        }
    }

    void add_with_multiply(unsigned i, unsigned j, Field lambda) {
        for (unsigned i1 = 0; i1 < N; i1++) {
            Field ch = lambda * _a[j][i1];
            _a[i][i1] += ch;
        }
    }

    void multiply_row(unsigned i, Field lambda) {
        for (unsigned i1 = 0; i1 < N; i1++) {
            _a[i][i1] *= lambda;
        }
    }

    Field Gauss(Matrix<M, N, Field> *Additional = nullptr) {
        Field mult_change = Field(1);
        unsigned i = 0, j = 0;
        while (i < M && j < N) {

            pair<Field, int> operation_result = _Gauss_operation(*this, i, j, Additional);
            int mov = operation_result.second;
            mult_change *= operation_result.first;
            if (mov == -1) { j++; }
            else {
                i++;
            }
        }
        return mult_change;
    }

    void Full_Blown_Gauss(Matrix<M, N, Field> *Additional = nullptr) {
        Gauss(Additional);
        unsigned i = 0, j = 0;
        while (i < M && j < N) {
            if (_a[i][j] == Field(0)) {
                j++;
            } else {
                Field st = _a[i][j];
                for (unsigned i1 = 0; i1 < i; i1++) {
                    Field lambda1 = Field(-1) * (_a[i1][j] / st);
                    add_with_multiply(i1, i, lambda1);
                    if (Additional != nullptr) {
                        Additional->add_with_multiply(i1, i, lambda1);
                    }
                }
                i++;
            }
        }
    }

    Field get_diagonal_mulltiply() const {
        Field mult = Field(1);
        for (unsigned i = 0; i < std::min(N, M); i++) {
            mult *= _a[i][i];
        }
        return mult;
    }

    Field det() const {
        Matrix<N, N, Field> A = *this;
        Field deter = A.Gauss();
        deter *= A.get_diagonal_mulltiply();
        return deter;
    }

    Matrix<N, M, Field> transposed() const {
        Matrix<N, M, Field> A;
        for (unsigned i = 0; i < N; i++) {
            for (unsigned j = 0; j < M; j++) {
                A[i][j] = _a[j][i];
            }
        }
        return A;
    }

    unsigned rank() const {
        Matrix<M, N, Field> A = *this;
        A.Gauss();
        unsigned rank = 0;
        unsigned i = 0, j = 0;
        while (i < M && j < N) {
            if (A._a[i][j] == Field(0)) {
                j++;
            } else {
                i++;
                rank++;
            }
        }
        return rank;
    }

    Matrix<N, N, Field> inverted() const {
        Matrix<N, N, Field> A = *this;
        A.invert();
        return A;
    };

    template<unsigned K>
    Matrix<M, K, Field> classic_multiply(const Matrix<N, K, Field> &B) const {
        Matrix<M, K, Field> ans;
        for (unsigned i = 0; i < M; i++) {
            for (unsigned j = 0; j < K; j++) {
                for (unsigned k = 0; k < N; k++) {
                    ans[i][j] += (_a[i][k] * B[k][j]);
                }
            }
        }
        return ans;
    }

    void invert() {
        My_Static_Assert<N == M>::assertion();
        Matrix<N, N, Field> A;
        for (unsigned i = 0; i < N; i++) {
            A._a[i][i] = Field(1);
        }
        Full_Blown_Gauss(&A);
        for (unsigned i = 0; i < N; i++) {
            A.multiply_row(i, Field(1) / _a[i][i]);
        }
        *this = A;
    };

    Field trace() const {
        My_Static_Assert<N == M>::assertion();
        Field tr = Field(0);
        for (unsigned i = 0; i < N; i++) {
            tr += _a[i][i];
        }
        return tr;
    }

    vector<Field> getRow(unsigned i) const {
        vector<Field> a1(N);
        for (unsigned j = 0; j < N; j++) {
            a1[j] = _a[i][j];
        }
        return a1;
    }

    vector<Field> getColumn(unsigned i) const {
        vector<Field> a1(N);
        for (unsigned j = 0; j < M; j++) {
            a1[j] = _a[j][i];
        }
        return a1;
    }

    template<unsigned K>
    Matrix<K, K, Field> _enlarge() const {
        Matrix<K, K, Field> A;
        for (unsigned i = 0; i < M; i++) {
            for (unsigned j = 0; j < N; j++) {
                A[i][j] = _a[i][j];
            }
        }
        return A;
    }

};

template<unsigned n, typename Field=Rational>
using SquareMatrix=Matrix<n, n, Field>;

template<unsigned N, typename Field=Rational>
const SquareMatrix<2 * N, Field> from_same_sized_blocks(const SquareMatrix<N, Field> a[2][2]) {
    SquareMatrix<2 * N, Field> M;
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++) {
            for (unsigned i1 = 0; i1 < 2; i1++) {
                for (unsigned j1 = 0; j1 < 2; j1++) {
                    M[i + i1 * N][j + j1 * N] = a[i1][j1][i][j];
                }
            }
        }
    }
    return M;
}

template<unsigned N, unsigned M, typename Field=Rational>
const Matrix<N, M, Field> operator+(const Matrix<N, M, Field> &a, const Matrix<N, M, Field> &b) {
    Matrix<N, M, Field> c = a;
    c += b;
    return c;
}

template<unsigned M, unsigned N, unsigned K, typename Field=Rational>
const Matrix<M, K, Field> strassen(Matrix<M, N, Field> &, Matrix<N, K, Field> &);

template<unsigned N, unsigned M, unsigned K, typename Field=Rational>
const Matrix<N, K, Field> operator*(const Matrix<N, M, Field> &a, const Matrix<M, K, Field> &b) {
    if (M <= 8 && N <= 8 && K <= 8) {
        return a.classic_multiply(b);
    }
    return strassen<N, M, K, Field>(a, b);
}

template<unsigned N, unsigned M, typename Field=Rational>
const Matrix<N, M, Field> operator-(const Matrix<N, M, Field> &a, const Matrix<N, M, Field> &b) {
    Matrix<N, M, Field> c = a;
    c -= b;
    return c;
}

template<unsigned N, unsigned M, typename Field=Rational>
const Matrix<N, M, Field> operator*(Field a, const Matrix<N, M, Field> &b) {
    Matrix<N, M, Field> c = b;
    c *= a;
    return c;
}

template<unsigned N, unsigned M, typename Field=Rational>
const Matrix<N, M, Field> operator*(const Matrix<N, M, Field> &b, Field a) {
    Matrix<N, M, Field> c = b;
    c *= a;
    return c;
}

template<unsigned N, unsigned M, typename Field=Rational>
bool operator==(const Matrix<N, M, Field> &a, const Matrix<N, M, Field> &b) {
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < M; j++) {
            if (a[i][j] != b[i][j]) {
                return false;
            }
        }
    }
    return true;
}

template<unsigned N, unsigned M, typename Field=Rational>
bool operator!=(const Matrix<N, M, Field> &a, const Matrix<N, M, Field> &b) {
    return !(a == b);
}

template<unsigned N, typename Field=Rational>
const SquareMatrix<N, Field> power_of_two_multiply(const SquareMatrix<N, Field> &A, const SquareMatrix<N, Field> &B) {
    if (N <= 16) {
        return A.classic_multiply(B);
    }
    My_Static_Assert<(N % 2 == 0) || N == 1>::assertion();
    SquareMatrix<N / 2, Field> ans[2][2];
    SquareMatrix<N / 2, Field> blocksA[2][2];
    SquareMatrix<N / 2, Field> blocksB[2][2];
    for (unsigned i = 0; i < N / 2; i++) {
        for (unsigned j = 0; j < N / 2; j++) {
            for (unsigned i1 = 0; i1 < 2; i1++) {
                for (unsigned j1 = 0; j1 < 2; j1++) {
                    blocksA[i1][j1][i][j] = A[i + i1 * N / 2][j + j1 * N / 2];
                    blocksB[i1][j1][i][j] = B[i + i1 * N / 2][j + j1 * N / 2];
                }
            }
        }
    }
    for (unsigned i = 0; i < 2; i++) {
        for (unsigned j = 0; j < 2; j++) {
            for (unsigned k = 0; k < 2; k++) {
                ans[i][j] += (power_of_two_multiply<N / 2, Field>(blocksA[i][k], blocksB[k][j]));
            }
        }
    }
    return from_same_sized_blocks<N / 2, Field>(ans);
}

template<unsigned N, typename Field=Rational>
const SquareMatrix<1, Field> power_of_two_multiply(SquareMatrix<1, Field> &A, SquareMatrix<1, Field> &B) {
    return A.classic_multiply(B);
}

template<unsigned M, unsigned N, unsigned K, typename Field>
const Matrix<M, K, Field> strassen(const Matrix<M, N, Field> &A, const Matrix<N, K, Field> &B) {
    static const unsigned two_degree = upper_degree_of_two_bound<static_max<static_max<M, N>::value, K>::value>::num;
    SquareMatrix<two_degree, Field> A1 = A.template _enlarge<two_degree>();
    SquareMatrix<two_degree, Field> B1 = B.template _enlarge<two_degree>();
    SquareMatrix<two_degree, Field> ANS = power_of_two_multiply<two_degree, Field>(A1, B1);
    Matrix<M, K, Field> ans;
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < K; j++) {
            ans[i][j] = ANS[i][j];
        }
    }
    return ans;
}



#endif //TEMPLATES_MATRIX_H

#endif
#endif //TEMPLATES_ZHOPA_H
