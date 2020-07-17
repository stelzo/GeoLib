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

#ifndef geo_h
#define geo_h

#include <string>
#include <vector>
#include <math.h>
#include <limits>

namespace geo
{
    // Converts degrees to radians.
    //
    // @param degrees
    // @return radians
    inline double deg2rad(double deg)
    {
        return deg * M_PI / 180;
    }

    // Converts radians to degrees.
    //
    // @param radians
    // @return degrees
    inline double rad2deg(double rad)
    {
        return rad * 180 / M_PI;
    }

    enum vectortype {
        rel,
        abs
    };
    // Vec2f represents a 2D vector in space which also can be a point.
    // It implements basic arithmetic and helper functions, to create and manipulate
    // points and vectors.
    //
    // For position information precision float was chosen because of its small memory footprint
    // and therefore better caching in algorithms and use cases where millimeter precision
    // is enough when 1 meter equals 1.0f.
    //
    // This class should never have any kind of setter methods. Algorithms and conversions need
    // to be implemented inside the library.
    class Vec2f
    {

    private:
        float _x;
        float _y;

        // Squared length for utilization in multiple algorithms
        //
        // @return length squared
        inline double _lengthSquared() const
        {
            return this->dot(*this);
        }

    public:
        // Constructor with initialization of both positions.
        //
        // @param x: position in x dimension
        // @param y: position in y dimension
        Vec2f(float x, float y);

        // Constructor initializes with x = y = 0.
        //
        Vec2f();

        // Destructor - does nothing in this case.
        //
        ~Vec2f();

        // Copy constructor initializes with x and y from v.
        //
        // @param vector to copy from
        Vec2f(const Vec2f &v);

        // Constructs a vector from a point to another point.
        //
        // @param from: origin
        // @param to: destination
        Vec2f(const Vec2f &from, const Vec2f &to);

        // Getter for x position.
        //
        // @return x position
        inline float x() const
        {
            return _x;
        }

        // Getter for y position.
        //
        // @return y position
        inline float y() const
        {
            return _y;
        }

        // Addition of both x and y positions.
        //
        // @param other vector
        // @return new vector with the values from the addition
        inline Vec2f operator+(const Vec2f &o) const
        {
            return Vec2f(_x + o._x, _y + o._y);
        }

        // Addition with another vector but the calling class gets the result.
        //
        // @param other vector
        // @return *this for concatenation
        inline Vec2f &operator+=(const Vec2f &o)
        {
            _x += o._x;
            _y += o._y;
            return *this;
        }

        // Assigns the calling class to the values of the given vector.
        //
        // @param other vector
        // @return *this for concatenation
        inline Vec2f &operator=(const Vec2f &o)
        {
            if (equals(o))
                return *this;
            _x = o._x;
            _y = o._y;
            return *this;
        }

        // Subtraction with another vector but the calling class gets the result.
        //
        // @param other vector
        // @return *this for concatenation
        inline Vec2f &operator-=(const Vec2f &o)
        {
            _x -= o._x;
            _y -= o._y;
            return *this;
        }

        // Subtraction with another vector.
        //
        // @param other vector
        // @return new vector with the result of the subtraction.
        inline Vec2f operator-(const Vec2f &o) const
        {
            return Vec2f(_x - o._x, _y - o._y);
        }

        // Unary inversion. All x and y change the sign.
        //
        // @return new vector with inverted direction
        inline Vec2f operator-() const
        {
            return Vec2f(-_x, -_y);
        }

        // Subtraction with a scalar. All elements get subtracted by s.
        // Equal to subtracting by a Vector [s, s].
        //
        // @param scalar
        inline Vec2f operator-(float s) const
        {
            return Vec2f(_x - s, _y - s);
        }

        // Division by a scalar. All elements get divided by s.
        //
        // @param scalar
        inline Vec2f operator/(float s) const
        {
            return s < std::numeric_limits<float>::epsilon() ? Vec2f() : Vec2f(_x / s, _y / s);
        }

        // Dot product with another vector. Mind the order. a.dot(b) == a . b.
        //
        // @param vector
        inline double dot(const Vec2f &v) const
        {
            return (_x * v._x) + (_y * v._y);
        }

        // Checks whether a vector has values near zero (std::limits::epsilon).
        //
        // @return true if values are near zero, else false
        inline bool zero() const
        {
            return fabs(_x) < std::numeric_limits<float>::epsilon() && fabs(_y) < std::numeric_limits<float>::epsilon();
        }

        // Calculates the projected point from the calling vector on vector v.
        // This is also the closest point possible to the destination on the line that v describes.
        //
        // @param vector to project on
        // @return point on v
        inline Vec2f projected_point(const Vec2f &v) const
        {
            return v.normalize() * dot(v.normalize());
        }

        // Calculates the direct vector to any point on v.
        // The resulting vector has the smallest length possible.
        //
        // @param vector that describes the line
        // @return vector from destination of calling vector to v
        inline Vec2f closest_vec2_to(const Vec2f &v) const
        {
            return Vec2f(*this, projected_point(v));
        }

        // Rotates a vector around a pivot (origin is default).
        // Positive radians mean clockwise,
        // negative radians mean anticlockwise rotation.
        //
        // @param radians to rotate
        // @param pivot to rotate around
        // @return rotated vector
        inline Vec2f rotate(float rad, const Vec2f &pivot = Vec2f(0, 0)) const
        {

            double s = sin(rad);
            double c = cos(rad);

            // translate back to origin and rotate clockwise
            float x = c * (_x - pivot._x) + (_y - pivot._y) * s;
            float y = s * -(_x - pivot._x) + (_y - pivot._y) * c;

            // translate back
            return Vec2f(x + pivot._x, y + pivot._y);
        }

        // Reflects a vector on another vector and returns the resulting bounce vector
        // in the normal direction of v.
        //
        // @param normal n
        // @return vector reflected on normal
        inline Vec2f reflect(const Vec2f &n) const
        {
            double f = this->dot(n);
            return Vec2f(_x - n._x * 2.0 * f, _y - n._y * 2.0 * f);
        }

        // Length of the vector in units.
        //
        // @return length
        inline double length() const
        {
            return sqrt(_lengthSquared());
        }

        // Normalizes vector of the calling vector and returns it.
        // Does not change the calling vector.
        //
        // @return normalized vector, length 1
        inline Vec2f normalize() const
        {
            double l = length();
            if (l < std::numeric_limits<float>::epsilon())
                return *this;
            return (*this / l);
        }

        // Scales the vector with a scalar s.
        //
        // @param scalar
        // @return scaled vector
        inline Vec2f operator*(float s) const
        {
            return Vec2f(_x * s, _y * s);
        }

        // Calculates the angle between the calling vector and v.
        // Positive results mean clockwise,
        // negative results mean anticlockwise rotation.
        //
        // Signing mechanism written by Patrick Hoffmann.
        //
        // @param v: vector/point to find angle to.
        // @param sign_coord: whether to check the x or y position for the
        // resulting sign.
        // @return angle in radians
        inline double angle_signed(const Vec2f &v, bool sign_coord_x = true) const
        {
            int sign = -1;
            if (sign_coord_x && v._x >= _x)
                sign = 1;
            if (!sign_coord_x && v._y >= _y)
                sign = 1;

            return sign * angle(v);
        }

        // Calculates the angle between the calling vector and v.
        //
        // @param vector to find angle to.
        // @return angle in radians
        inline double angle(const Vec2f &v) const
        {
            double len_a = length();
            double len_b = v.length();

            if (len_a < 0.00001 || len_b < 0.00001)
                return 0.0;

            double angle = acos(dot(v) / (len_a * len_b));

            return angle < 0.00001 ? 0.0f : angle;
        }

        // Check equality of two vectors. Same functionality as operator==
        // but mainly used for testing.
        //
        // @param vector v to check equality against.
        // @return true if same values, else false
        inline bool equals(const Vec2f &v) const
        {
            return fabs(_x - v._x) < 0.00001 && fabs(_y - v._y) < 0.00001;
        }

        // Check equality of two vectors.
        //
        // @param vector v to check equality against.
        // @return true if same values, else false
        inline bool operator==(const Vec2f &b) const
        {
            return equals(b);
        }

        // Check inequality of two vectors.
        //
        // @param vector v to check equality against.
        // @return true if different values, else false
        inline bool operator!=(const Vec2f &b) const
        {
            return !equals(b);
        }

        // Implicit string conversion. Minimal string representation
        // of a vector. Same functionality as to_string().
        //
        // @return string representation of the vector.
        operator std::string() const;

        // Minimal string representation of a vector.
        //
        // @return string representation of the vector.
        std::string to_string() const;
    };

    // Multiple points can be connected to define an object.
    // The main structure describing that object is called a polygon.
    //
    // Landau notation for all functions is given.
    // Vector intern reallocation is ignored.
    // n is scaling linearly with the amount of points in the polygon
    // m is scaling linearly with the total amount of points passed
    //
    // !!! A polygon only knows points, so each Vec2f will be interpreted as a point. !!!
    // Points will be called vertex or vertices.
    //
    // The class manages the points and gives access to general manipulations of the structure.
    // For specific functions, derive from this class.
    class Polygon2
    {
    private:
        std::vector<Vec2f> vertices;
        vectortype type = abs;

    public:
        Polygon2() = default;

        // Builds a polygon with the possibility to accept relative vectors.
        // O(1) time.
        //
        // @param polygontype what type of vectors this polygon should accept
        Polygon2(vectortype polygontype);

        // Builds a polygon from all sides of a thing by combining the sides.
        // O(m) time.
        // 
        // @param sides of a structure described in points
        Polygon2(const std::vector<std::vector<Vec2f>> &sides);

        // Builds a polygon from all sides of a thing by combining the sides with the possibility to accept relative vectors.
        // O(m) time.
        //
        // @param sides of a structure described in points
        // @param polygontype what type of vectors this polygon should accept
        Polygon2(const std::vector<std::vector<Vec2f>> &sides, vectortype polygontype);

        // Builds a polygon from points.
        // O(m) time.
        //
        // @param vector of all points.
        Polygon2(const std::vector<Vec2f> &vertices);

        // Builds a polygon from points with the possibility to accept relative vectors.
        // O(m) time.
        //
        // @param vector of all points.
        // @param polygontype what type of vectors this polygon should accept
        Polygon2(const std::vector<Vec2f> &vertices, vectortype polygontype);

        // Add a single vertex to the vertices programmatically.
        // O(1) time.
        //
        // @param vertex to add to the vertices.
        inline void add_vertex(const Vec2f &vertex)
        {
            switch (type) {
                case rel:
                    vertices.push_back(vertices.back() + vertex);
                    break;
                case abs:
                    vertices.push_back(vertex);
                    break;
            }
        }

        // Checks whether a point is inside a polygon.
        // O(n) time.
        //
        // @param point to check
        // @return true if inside, else false
        bool contains(const Vec2f &p);

        // Smoothes the polygon with catmull-rom splines.
        // [WARNING] This function is not tested (because it is hard to do).
        // O(n*(1/distance)) time.
        //
        // @param alpha type of spline (0.0 = uniform, 0.5 = centripetal, 1.0 = chordal). Must be between 0.0 and 1.0
        // @param tension how tight the spline should be, 1.0 would result in linear splines. Must be between 0.0 and 1.0.
        // @param distance distance between new points in 1/100 of the distance between to original points
        void smooth(float alpha, float tension, unsigned int distance);

        // Adds a vertex to the polygon.
        // O(n) time complexity, because it returns a copy of the polygon.
        //
        // @param vertex to add
        // @return new polygon with added vertex
        inline Polygon2 operator+(const Vec2f &vertex) const
        {
            Polygon2 ret;
            switch (type) {
                case rel:
                    ret = Polygon2({vertices, {vertices.back() + vertex}});
                    break;
                case abs:
                    ret = Polygon2({vertices, {vertex}});
                    break;
            }
            return ret;
        }

        // Adds a vertex to the polygon but inclusive and pretty.
        // O(1) time.
        //
        // @param vertex to add
        // @return new polygon with added vertex
        inline Polygon2 &operator<<(const Vec2f &vertex)
        {
            return (*this += vertex);
        }

        // Adds a vertex to the polygon but inclusive.
        // O(1) time.
        //
        // @param vertex to add
        // @return *this for concatenation
        inline Polygon2 &operator+=(const Vec2f &vertex)
        {
            switch (type) {
                case rel:
                    add_vertex(vertices.back() + vertex);
                    break;
                case abs:
                    add_vertex(vertex);
                    break;
            }
            return *this;
        }

        // Moves all points of the Polygon.
        // O(n) time.
        //
        // @param movement of the polygon
        void move(const Vec2f& movement);

        // Rotates all points of the Polygon.
        // Positive radians mean clockwise,
        // negative radians mean anticlockwise rotation.
        // O(n) time.
        //
        // @param radians to rotate
        // @param pivot to rotate around
        void rotate(float rad, const Vec2f& pivot = Vec2f());

        inline void setType(vectortype polygontype)
        {
            type = polygontype;
        }

        // Sizes the polygon along normals to allow for deadzones
        // Positive values mean enlargement.
        // Negative values mean shrinking.
        // O(n²)
        //
        // @param dist how far the polygon should be scaled
        Polygon2 size(float dist);


        // Minimal string representation of a polygon.
        // O(n)
        //
        // @return string representation of the polygon.
        std::string to_string() const;

        // Copying a polygon
        // O(n)
        //
        // @param clone the polygon to clone
        Polygon2 &operator=(const Polygon2& clone) {
            vertices = clone.vertices;
            type = clone.type;
            return *this;
        }

        // Copy constructor.
        // O(m) time.
        //
        // @param polygon to copy from.
        Polygon2(const Polygon2 &v);
    };

} // namespace geo

#endif // geo_h