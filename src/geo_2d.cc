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

Vec2f::~Vec2f() {}

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

Polygon2::Polygon2(const std::vector<std::vector<Vec2f>> &sides)
{
    for (auto side : sides)
    {
        vertices.insert(vertices.end(), side.begin(), side.end());
    }
}

Polygon2::Polygon2(const std::vector<Vec2f> &vertices) : vertices(vertices) {}

void Polygon2::add_vertex(const Vec2f &vertex)
{
    vertices.push_back(vertex);
}

bool Polygon2::is_inside(const Vec2f &p)
{
    size_t nvert = vertices.size();
    size_t i, j;
    i = j = 0;
    bool c = false;

    for (i = 0, j = nvert - 1; i < nvert; j = i++)
    {
        if (((vertices[i].y() > p.y()) != (vertices[j].y() > p.y())) &&
            (p.x() < (vertices[j].x() - vertices[i].x()) * (p.y() - vertices[i].y()) / (vertices[j].y() - vertices[i].y()) + vertices[i].x()))
        {
            c = !c;
        }
    }

    return c;
}

void Polygon2::smooth(float alpha, float tension, unsigned int distance)
{
    size_t size = vertices.size();
    std::vector<Vec2f> smoothed;

    for (unsigned int i = 0; i < size - 3; i++)
    {
        float t01 = pow(Vec2f(vertices[i], vertices[i + 1]).length(), alpha);
        float t12 = pow(Vec2f(vertices[i + 1], vertices[i + 2]).length(), alpha);
        float t23 = pow(Vec2f(vertices[i + 2], vertices[i + 3]).length(), alpha);

        float x1, y1, x2, y2;

        x1 = (1.0 - tension) * (vertices[i + 2].x() - vertices[i + 1].x() + t12 * ((vertices[i + 1].x() - vertices[i].x()) / t01 - (vertices[i + 2].x() - vertices[i].x()) / (t01 + t12)));
        y1 = (1.0 - tension) * (vertices[i + 2].y() - vertices[i + 1].y() + t12 * ((vertices[i + 1].y() - vertices[i].y()) / t01 - (vertices[i + 2].y() - vertices[i].y()) / (t01 + t12)));

        x2 = (1.0 - tension) * (vertices[i + 2].x() - vertices[i + 1].x() + t12 * ((vertices[i + 3].x() - vertices[i + 2].x()) / t23 - (vertices[i + 3].x() - vertices[i + 1].x()) / (t12 + t23)));
        y2 = (1.0 - tension) * (vertices[i + 2].y() - vertices[i + 1].y() + t12 * ((vertices[i + 3].y() - vertices[i + 2].y()) / t23 - (vertices[i + 3].y() - vertices[i + 1].y()) / (t12 + t23)));

        float _a_x, _a_y, _b_x, _b_y, _c_x, _c_y, _d_x, _d_y;
        _a_x = 2.0 * (vertices[i + 1].x() - vertices[i + 2].x()) + x1 + x2;
        _a_y = 2.0 * (vertices[i + 1].y() - vertices[i + 2].y()) + y1 + y2;
        _b_x = -3.0 * (vertices[i + 1].x() - vertices[i + 2].x()) - x1 - x1 - x2;
        _b_y = -3.0 * (vertices[i + 1].y() - vertices[i + 2].y()) - y1 - y1 - y2;
        _c_x = x1;
        _c_y = y1;
        _d_x = vertices[i + 1].x();
        _d_y = vertices[i + 1].y();

        for (int j = 0; j < 100; j = j + distance)
        {
            float t = (float)j / 100;
            float x_s = _a_x * t * t * t + _b_x * t * t + _c_x * t + _d_x;
            float y_s = _a_y * t * t * t + _b_y * t * t + _c_y * t + _d_y;

            smoothed.emplace_back(x_s, y_s);
        }
    }

    vertices = smoothed;
}

Polygon2 Polygon2::operator+(const Vec2f &vertex) const
{
    return Polygon2({vertices, {vertex}});
}

Polygon2 &Polygon2::operator+=(const Vec2f &vertex)
{
    add_vertex(vertex);
    return *this;
}

Polygon2 &Polygon2::operator<<(const Vec2f &vertex)
{
    return (*this += vertex);
}

Polygon2::Polygon2(const Polygon2 &v) : vertices(v.vertices) {}