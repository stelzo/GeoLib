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

Vec2f::Vec2f(const Vec2f &from, const Vec2f &to) : _x(to._x - from._x), _y(to._y - from._y) {}

Vec2f Vec2f::projected_point(const Vec2f &v)
{
    return v.normalize() * dot(v.normalize());
}

Vec2f Vec2f::closest_vec2_to(const Vec2f &v)
{
    return Vec2f(*this, projected_point(v));
}

bool Vec2f::equals(const Vec2f &v) const
{
    return fabs(_x - v._x) < 0.00001 && fabs(_y - v._y) < 0.00001;
}

bool Vec2f::operator==(const Vec2f &b) const
{
    return equals(b);
}

std::string Vec2f::to_string() const
{
    std::stringstream ss;
    ss << "[Vec2] "
       << "x: " << _x << " -- "
       << "y: " << _y;

    return ss.str();
}

Vec2f::operator std::string() const
{
    return to_string();
}

bool Vec2f::operator!=(const Vec2f &b) const
{
    return !equals(b);
}

bool Vec2f::empty() const
{
    return _x < std::numeric_limits<float>::epsilon() && _y < std::numeric_limits<float>::epsilon();
}

Vec2f Vec2f::rotate(float rad, const Vec2f &pivot)
{

    double s = sin(rad);
    double c = cos(rad);

    // translate back to origin and rotate clockwise
    float x = c * (_x - pivot._x) + (_y - pivot._y) * s;
    float y = s * -(_x - pivot._x) + (_y - pivot._y) * c;

    // translate back
    return Vec2f(x + pivot._x, y + pivot._y);
}

Vec2f Vec2f::reflect(const Vec2f &v, const Vec2f &n)
{
    float f = _x * n._x + _y * n._y;
    return Vec2f(_x - n._x * 2.0 * f, _y - n._y * 2.0 * f);
}

Vec2f::Vec2f(float x, float y) : _x(x), _y(y) {}

Vec2f::Vec2f() : _x(0), _y(0) {}

Vec2f::Vec2f(const Vec2f &v) : _x(v._x), _y(v._y) {}

float Vec2f::x() const
{
    return _x;
}

float Vec2f::y() const
{
    return _y;
}

Vec2f Vec2f::operator+(const Vec2f &o) const
{
    return Vec2f(_x + o._x, _y + o._y);
}

Vec2f &Vec2f::operator+=(const Vec2f &o)
{
    _x += o._x;
    _y += o._y;
    return *this;
}

Vec2f &Vec2f::operator=(const Vec2f &o)
{
    if (equals(o))
        return *this;
    _x = o._x;
    _y = o._y;
    return *this;
}

Vec2f Vec2f::operator-()
{
    return Vec2f(-_x, -_y);
}

Vec2f &Vec2f::operator-=(const Vec2f &o)
{
    _x -= o._x;
    _y -= o._y;
    return *this;
}

Vec2f Vec2f::operator-(const Vec2f &o) const
{
    return Vec2f(_x - o._x, _y - o._y);
}

Vec2f Vec2f::operator-(float s) const
{
    return Vec2f(_x - s, _y - s);
}

Vec2f Vec2f::operator/(float s) const
{
    return s < std::numeric_limits<float>::epsilon() ? Vec2f() : Vec2f(_x / s, _y / s);
}

double Vec2f::dot(const Vec2f &v) const
{
    return (_x * v._x) + (_y * v._y);
}

double Vec2f::lengthSquared() const
{
    return this->dot(*this);
}

double Vec2f::length() const
{
    return sqrt(lengthSquared());
}

Vec2f Vec2f::normalize() const
{
    double l = length();
    if (l < std::numeric_limits<float>::epsilon())
        return *this;
    return (*this / l);
}

Vec2f Vec2f::operator*(float s) const
{
    return Vec2f(_x * s, _y * s);
}

// Based on a function written by Patrick Hoffmann
double Vec2f::angle(const Vec2f &v, bool sign_by_x_coord)
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
