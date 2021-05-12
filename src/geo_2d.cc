#include <stdio.h>
#include "geo.h"
#include <math.h>

#include <sstream>

using namespace geo;

Vec2f::Vec2f(const Vec2f &from, const Vec2f &to) : _x(to._x - from._x), _y(to._y - from._y) {}

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

Vec2f::Vec2f(float x, float y) : _x(x), _y(y) {}

Vec2f::Vec2f() : _x(0), _y(0) {}

Vec2f::Vec2f(const Vec2f &v) : _x(v._x), _y(v._y) {}

Polygon2::Polygon2(const std::vector<std::vector<Vec2f>> &sides)
{
    for (auto side : sides)
    {
        vertices.insert(vertices.end(), side.begin(), side.end());
    }
}

Polygon2::Polygon2(const std::vector<Vec2f> &vertices) : vertices(vertices) {}

bool Polygon2::contains(const Vec2f &p)
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

void Polygon2::move(const Vec2f& movement)
{
    for (auto& vertex : vertices) {
        vertex += movement;
    }
}

void Polygon2::rotate(float rad, const Vec2f& pivot)
{
    for (auto& vertex : vertices) {
        vertex = vertex.rotate(rad, pivot);
    }
}

Polygon2::Polygon2(const Polygon2 &v) : vertices(v.vertices) {}

Polygon2::Polygon2 (vectortype polygontype) {
    type = polygontype;
}

Polygon2::Polygon2 (const std::vector<std::vector<Vec2f>> &sides, vectortype polygontype) {
    for (auto side : sides)
    {
        vertices.insert(vertices.end(), side.begin(), side.end());
    }
    type = polygontype;
}

Polygon2::Polygon2 (const std::vector<Vec2f> &vertices, vectortype polygontype) : vertices(vertices) {
    type = polygontype;
}

Polygon2 Polygon2::sized (float dist) {
    std::vector<Vec2f> sized;

    Vec2f con(vertices.back(), vertices[1]);
    Vec2f rn(con.y(), -con.x());
    sized.push_back(rn.normalize() * dist + vertices[0]);

    for (int i = 1; i < (int) vertices.size() - 1; ++i) {
        Vec2f connection(vertices[i - 1], vertices[i + 1]);
        Vec2f rnorm(connection.y(), -connection.x());
        sized.push_back(rnorm.normalize() * dist + vertices[i]);
    }


    Vec2f c(vertices[vertices.size() - 2], vertices[0]);
    Vec2f r(c.y(), -c.x());
    sized.push_back(r.normalize() * dist + vertices[vertices.size() - 1]);

    return Polygon2(sized);
}

Polygon2 Polygon2::scaled (float dist) {
    std::vector<Vec2f> sized;

    Vec2f con(vertices.back(), vertices[1]);
    Vec2f rn(con.y(), -con.x());
    sized.push_back(rn.normalize() * dist + vertices[0]);

    for (int i = 1; i < (int) vertices.size() - 1; ++i) {
        Vec2f connection(vertices[i - 1], vertices[i + 1]);
        Vec2f rnorm(connection.y(), -connection.x());
        sized.push_back(rnorm.normalize() * dist + vertices[i]);
    }


    Vec2f c(vertices[vertices.size() - 2], vertices[0]);
    Vec2f r(c.y(), -c.x());
    sized.push_back(r.normalize() * dist + vertices[vertices.size() - 1]);

    return Polygon2(sized);
}

std::string Polygon2::to_string () const {
    std::stringstream ss;
    ss << "[Polygon2]\n";
    for (auto vertex : vertices) {
        ss << vertex.to_string() << "\n";
    }

    return ss.str();
}

double Polygon2::area () {
    double phalf = 0, nhalf = 0;

    for (int i = 0; i < (int) vertices.size(); ++i) {
        phalf += vertices[i].x() * vertices[(i + 1) % vertices.size()].y();
        nhalf += vertices[(i + 1) % vertices.size()].x() * vertices[i].y();
    }

    return std::abs(phalf - nhalf) / 2;
}
