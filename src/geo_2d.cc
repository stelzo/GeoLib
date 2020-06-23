//
//  geo_2d.cpp
//  libgeo
//
//  Created by Christopher Sieh on 21.06.20.
//  Copyright © 2020 steado. All rights reserved.
//

#include <stdio.h>
#include "geo_2d.h"
#include <math.h>
#include <limits>

#include <sstream>

using namespace geo;

template <class _T>
Vec2<_T>::Vec2(const Vec2<_T> &from, const Vec2<_T> &to) : _x(to._x - from._x), _y(to._y - from._x) {}

template <class _T>
Vec2<_T> Vec2<_T>::projected_point(const Vec2<_T> &v)
{
    return v.normalize() * dot(v.normalize());
}

template <class _T>
Vec2<_T> Vec2<_T>::closest_vec2_to(const Vec2<_T> &v)
{
    return Vec2<_T>(*this, projected_point(v));
}

// TODO std::epsilon is too small
template <class _T>
bool Vec2<_T>::equals(const Vec2<_T> &v) const
{
    return std::abs(_x-v._x) < 0.00001 && std::abs(_y-v._y) < 0.00001;
}

template <class _T>
bool Vec2<_T>::operator==(const Vec2<_T> &b) const
{
    return equals(b);
}

template <class _T>
std::string Vec2<_T>::to_string() const
{
    std::stringstream ss;
    ss << "[Vec2] " << "x: " << _x << " -- " << "y: " << _y;

    return ss.str();
}

template <class _T>
Vec2<_T>::operator std::string() const
{
    return to_string();
}

template <class _T>
bool Vec2<_T>::operator!=(const Vec2<_T> &b) const
{
    return !equals(b);
}

template <class _T>
bool Vec2<_T>::empty() const
{
    return _x < std::numeric_limits<_T>::epsilon() && _y < std::numeric_limits<_T>::epsilon();
}

template <class _T>
Vec2<_T> Vec2<_T>::rotate(float rad, const Vec2<_T> &pivot)
{

    double s = sin(rad);
    double c = cos(rad);

    // translate back to origin and rotate clockwise
    float x = c * (_x-pivot._x) + (_y-pivot._y) * s;
    float y = s * -(_x-pivot._x) + (_y-pivot._y) * c;

    // translate back
    return Vec2<_T>(x + pivot._x, y + pivot._y);
}

template <class _T>
Vec2<_T> Vec2<_T>::reflect(const Vec2<_T> &v, const Vec2<_T> &n)
{
    float f = _x * n._x + _y * n._y;
    return Vec2<_T>(_x - n._x * 2.0 * f, _y - n._y * 2.0 * f);
}

template <class _T>
Vec2<_T>::Vec2(_T x, _T y) : _x(x), _y(y) {}

template <class _T>
Vec2<_T>::Vec2() : _x(0), _y(0) {}

template <class _T>
Vec2<_T>::Vec2(const Vec2<_T> &v) : _x(v._x), _y(v._y) {}

template <class _T>
_T Vec2<_T>::x() const
{
    return _x;
}

template <class _T>
_T Vec2<_T>::y() const
{
    return _y;
}

template <class _T>
Vec2<_T> Vec2<_T>::operator+(const Vec2<_T> &o) const
{
    return Vec2<_T>(_x + o._x, _y + o._y);
}

template <class _T>
Vec2<_T> &Vec2<_T>::operator+=(const Vec2<_T> &o)
{
    _x += o._x;
    _y += o._y;
    return *this;
}

template <class _T>
Vec2<_T> &Vec2<_T>::operator=(const Vec2<_T> &o)
{
    if (equals(o))
        return *this;
    _x = o._x;
    _y = o._y;
    return *this;
}

template <class _T>
Vec2<_T> Vec2<_T>::operator-()
{
    return Vec2(-_x, -_y);
}

template <class _T>
Vec2<_T> &Vec2<_T>::operator-=(const Vec2<_T> &o)
{
    _x -= o._x;
    _y -= o._y;
    return *this;
}

template <class _T>
Vec2<_T> Vec2<_T>::operator-(const Vec2<_T> &o) const
{
    return Vec2<_T>(_x - o._x, _y - o._y);
}

template <class _T>
Vec2<_T> Vec2<_T>::operator-(_T s) const
{
    return Vec2<_T>(_x - s, _y - s);
}

template <class _T>
Vec2<_T> Vec2<_T>::operator/(_T s) const
{
    return s < std::numeric_limits<_T>::epsilon() ? Vec2<_T>() : Vec2<_T>(_x / s, _y / s);
}

template <class _T>
double Vec2<_T>::dot(const Vec2<_T> &v) const
{
    return (_x * v._x) + (_y * v._y);
}

template <class _T>
double Vec2<_T>::lengthSquared() const
{
    return this->dot(*this);
}

template <class _T>
double Vec2<_T>::length() const
{
    return sqrt(lengthSquared());
}

template <class _T>
Vec2<_T> Vec2<_T>::normalize() const
{
    double l = length();
    if (l < std::numeric_limits<_T>::epsilon())
        return *this;
    return (*this / l);
}

template <class _T>
Vec2<_T> Vec2<_T>::operator*(_T s) const
{
    return Vec2<_T>(_x * s, _y * s);
}

// Based on a function written by Patrick Hoffmann
template <class _T>
double Vec2<_T>::angle(const Vec2<_T> &v, bool sign_by_x_coord)
{
    int sign = -1;
    if (sign_by_x_coord && v._x >= _x)
        sign = 1;
    if (!sign_by_x_coord && v._y >= _y)
        sign = 1;

    double len_a = length();
    double len_b = v.length();

    if (len_a < 0.00001 || len_b < 0.00001)
        return 0.0;

    double angle = sign * acos(this->dot(v) / (len_a * len_b));

    return angle < 0.00001 ? 0.0f : angle;
}
