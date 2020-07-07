// GeoLib - 2D/3D Geometry
// Copyright (C) 2020 Christopher Sieh (stelzo@steado.de)

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include "geo.h"
#include <math.h>
#include <limits>

#include <sstream>

using namespace geo;

Vec2f::Vec2f(const Vec2f &from, const Vec2f &to) : _x(to._x - from._x), _y(to._y - from._y) {}

Vec2f Vec2f::projected_point(const Vec2f &v) const
{
    return v.normalize() * dot(v.normalize());
}

Vec2f Vec2f::closest_vec2_to(const Vec2f &v) const
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

bool Vec2f::zero() const
{
    return fabs(_x) < std::numeric_limits<float>::epsilon() && fabs(_y) < std::numeric_limits<float>::epsilon();
}

Vec2f Vec2f::rotate(float rad, const Vec2f &pivot) const
{

    double s = sin(rad);
    double c = cos(rad);

    // translate back to origin and rotate clockwise
    float x = c * (_x - pivot._x) + (_y - pivot._y) * s;
    float y = s * -(_x - pivot._x) + (_y - pivot._y) * c;

    // translate back
    return Vec2f(x + pivot._x, y + pivot._y);
}

Vec2f Vec2f::reflect(const Vec2f &n) const
{
    double f = this->dot(n);
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

Vec2f Vec2f::operator-() const
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

double Vec2f::_lengthSquared() const
{
    return this->dot(*this);
}

double Vec2f::length() const
{
    return sqrt(_lengthSquared());
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

double Vec2f::angle_signed(const Vec2f &v, bool sign_coord_x) const
{
    int sign = -1;
    if (sign_coord_x && v._x >= _x)
        sign = 1;
    if (!sign_coord_x && v._y >= _y)
        sign = 1;

    return sign * angle(v);
}

double Vec2f::angle(const Vec2f &v) const
{
    double len_a = length();
    double len_b = v.length();

    if (len_a < 0.00001 || len_b < 0.00001)
        return 0.0;

    double angle = acos(dot(v) / (len_a * len_b));

    return angle < 0.00001 ? 0.0f : angle;
}
