// GeoLib - 2D/3D Geometry
// Copyright (c) 2020-present Christopher Sieh and other contributors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef geo_h
#define geo_h

#include <string>
#include <utility>
#include <vector>
#include <math.h>
#include <limits>
#include <algorithm>

namespace geo {
// Converts degrees to radians.
//
// @param degrees
// @return radians
inline double deg2rad (double deg) {
    return deg * M_PI / 180;
}

// Converts radians to degrees.
//
// @param radians
// @return degrees
inline double rad2deg (double rad) {
    return rad * 180 / M_PI;
}

enum vectortype {
    abs = 0,
    rel = 1
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
class Vec2f {

protected:
    float _x;
    float _y;

private:
    // Squared length for utilization in multiple algorithms
    //
    // @return length squared
    inline double _lengthSquared () const {
        return this->dot(*this);
    }

public:
    // Constructor with initialization of both positions.
    //
    // @param x: position in x dimension
    // @param y: position in y dimension
    Vec2f (float x, float y);

    // Constructor initializes with x = y = 0.
    //
    Vec2f ();

    // Destructor - does nothing in this case.
    //
    ~Vec2f ();

    // Copy constructor initializes with x and y from v.
    //
    // @param vector to copy from
    Vec2f (const Vec2f &v);

    // Constructs a vector from a point to another point.
    //
    // @param from: origin
    // @param to: destination
    Vec2f (const Vec2f &from, const Vec2f &to);

    // Getter for x position.
    //
    // @return x position
    inline float x () const {
        return _x;
    }

    // Getter for y position.
    //
    // @return y position
    inline float y () const {
        return _y;
    }

    // Addition of both x and y positions.
    //
    // @param other vector
    // @return new vector with the values from the addition
    inline Vec2f operator+ (const Vec2f &o) const {
        return Vec2f(_x + o._x, _y + o._y);
    }

    // Addition with another vector but the calling class gets the result.
    //
    // @param other vector
    // @return *this for concatenation
    inline Vec2f &operator+= (const Vec2f &o) {
        _x += o._x;
        _y += o._y;
        return *this;
    }

    // Assigns the calling class to the values of the given vector.
    //
    // @param other vector
    // @return *this for concatenation
    inline Vec2f &operator= (const Vec2f &o) {
        if (equals(o)) {
            return *this;
        }
        _x = o._x;
        _y = o._y;
        return *this;
    }

    // Subtraction with another vector but the calling class gets the result.
    //
    // @param other vector
    // @return *this for concatenation
    inline Vec2f &operator-= (const Vec2f &o) {
        _x -= o._x;
        _y -= o._y;
        return *this;
    }

    // Subtraction with another vector.
    //
    // @param other vector
    // @return new vector with the result of the subtraction.
    inline Vec2f operator- (const Vec2f &o) const {
        return Vec2f(_x - o._x, _y - o._y);
    }

    // Unary inversion. All x and y change the sign.
    //
    // @return new vector with inverted direction
    inline Vec2f operator- () const {
        return Vec2f(-_x, -_y);
    }

    // Subtraction with a scalar. All elements get subtracted by s.
    // Equal to subtracting by a Vector [s, s].
    //
    // @param scalar
    inline Vec2f operator- (float s) const {
        return Vec2f(_x - s, _y - s);
    }

    // Division by a scalar. All elements get divided by s.
    //
    // @param scalar
    inline Vec2f operator/ (float s) const {
        return s < std::numeric_limits<float>::epsilon() ? Vec2f() : Vec2f(_x / s, _y / s);
    }

    // Dot product with another vector. Mind the order. a.dot(b) == a . b.
    //
    // @param vector
    inline double dot (const Vec2f &v) const {
        return (_x * v._x) + (_y * v._y);
    }

    // Checks whether a vector has values near zero (std::limits::epsilon).
    //
    // @return true if values are near zero, else false
    inline bool zero () const {
        return fabs(_x) < std::numeric_limits<float>::epsilon() &&
               fabs(_y) < std::numeric_limits<float>::epsilon();
    }

    // Calculates the projected point from the calling vector on vector v.
    // This is also the closest point possible to the destination on the line that v describes.
    //
    // @param vector to project on
    // @return point on v
    inline Vec2f projected_point (const Vec2f &v) const {
        return v.normalize() * dot(v.normalize());
    }

    // Calculates the direct vector to any point on v.
    // The resulting vector has the smallest length possible.
    //
    // @param vector that describes the line
    // @return vector from destination of calling vector to v
    inline Vec2f closest_vec2_to (const Vec2f &v) const {
        return Vec2f(*this, projected_point(v));
    }

    // Calculates the closest vector on v.
    // The resulting vector has the smallest distance to any point on v.
    // [WARNING] This function is not tested.
    //
    // @param vector that describes the line
    // @return vector to a point on v that is closest this vector
    inline Vec2f closest_connection (const Vec2f &v) const {
        Vec2f  proj   = projected_point(v);
        Vec2f  x(proj.x(), v.x()), y(proj.y(), v.y());
        double x_comp = proj.x() / v.x(), y_comp = proj.y() / v.y();

        // cant use in() because of three output cases
        if ((y_comp < 0 && x_comp < 0) || (x.zero() && y_comp < 0) || (y.zero() && x_comp < 0)) {
            return -*this; // to [0, 0]
        }
        else if ((y_comp < 1 && x_comp < 1) || (x.zero() && y_comp < 1) || (y.zero() && x_comp < 1)) {
            return Vec2f(*this, proj);
        }
        else {
            return Vec2f(*this, v);
        }

    }

    // Finds the closest connection to a Spline between Points.
    // [WARNING] The given vector must be made of Points, not a Spline.
    //
    // @param points the Points on which Spline the closest connection should be found
    // @return the closest point on the Spline
    inline Vec2f closest_to (std::vector<Vec2f> points) {

        // create Spline and set up defaults
        std::vector<Vec2f> spline        = to_spline(points);
        double             smallest_dist = MAXFLOAT;
        Vec2f              closest_point;

        // Iterate through all possible Segments
        for (unsigned int j = 0; j < points.size() - 1; j++) {
            Vec2f compare(points[j], *this); //possible segfault?

            Vec2f candidate = compare.closest_connection(spline[j]);

            double p_proj_d = candidate.length();

            // Find the closest point
            if (p_proj_d < smallest_dist) {
                closest_point = candidate;
                smallest_dist = p_proj_d;
            }
        }

        return closest_point;
    }

    // Finds the closest point on a Spline between Points.
    // [WARNING] The given vector must be made of Points, not a Spline.
    //
    // @param points the Points on which Spline the closest point should be found
    // @return the closest point on the Spline
    inline Vec2f closest_on (std::vector<Vec2f> points) {

        return *this + closest_to(std::move(points));
    }

    // Creates a Spline from a given set of points
    // [WARNING] This function is only indirectly tested.
    //
    // @param points the points to be made into a spline
    // @return the created spline
    friend std::vector<Vec2f> to_spline (std::vector<Vec2f> points) {
        std::vector<Vec2f> spline;
        for (unsigned int  i = 0; i < points.size() - 1; i++) {
            spline.emplace_back(points[i], points[i + 1]);
        }
        return spline;
    }

    // Finds the intersection between this vector and a direction from an origin.
    //
    // @param origin the Point where dir originates
    // @param dir the normalized vector in the desired direction
    // @return the point where the vectors intersect, [0, 0] if parallel
    inline Vec2f intersection (const Vec2f& origin, const Vec2f& dir) const {

        if (this->normalize().dot(dir) == 1) {
            // Parallel
            return Vec2f(0, 0);
        }

        Vec2f con  = origin.closest_vec2_to(*this);
        Vec2f proj = dir.projected_point(con);


        // use dot product only for determining opposing vectors
        double scale = con.dot(proj) * con.length() / proj.length();
        //printf("Scale: %f, Direction: %s, Dot: %f\n", scale, (dir * scale).to_string().c_str(), con.dot(proj));
        return origin + (dir * scale);
    }

    // Finds the intersection between this vector on a direction from an origin.
    //
    // @param origin the Point where dir originates
    // @param dir the vector in the desired direction
    // @return the point where the vectors intersect, or [0, 0] if not intersecting
    inline Vec2f intersection_on (const Vec2f& origin, const Vec2f& dir) const {

        // get the intersection point
        Vec2f intersection = this->intersection(origin, dir.normalize());

        // check if on dir
        if (Vec2f(origin, intersection).in(dir)) {
            return intersection;
        }
        return Vec2f(0, 0);
    }


    // Finds the intersection on a Spline between Points.
    // Casts a ray onto a Spline between Points.
    // [WARNING] The given vector must be made of Points, not a Spline.
    //
    // @param points the Points on which Spline the intersection should be found
    // @return the intersection with the Spline, or [0, 0] if not intersecting
    inline Vec2f intersection_on (std::vector<Vec2f> points) const {

        // create Spline and set up defaults
        std::vector<Vec2f> spline        = to_spline(points);
        double             smallest_dist = MAXFLOAT;
        Vec2f              closest_point;

        // Iterate through all possible Segments
        for (unsigned int j = 0; j < points.size() - 1; j++) {

            Vec2f candidate = this->intersection_on(points[j], spline[j]);

            if (!candidate.zero()) {
                double p_proj_d = candidate.length();

                // Find the closest point
                if (p_proj_d < smallest_dist) {
                    closest_point = candidate;
                    smallest_dist = p_proj_d;
                }
            }
        }

        return closest_point;
    }

    // Rotates a vector around a pivot (origin is default).
    // Positive radians mean clockwise,
    // negative radians mean anticlockwise rotation.
    //
    // @param radians to rotate
    // @param pivot to rotate around
    // @return rotated vector
    inline Vec2f rotate (float rad, const Vec2f &pivot = Vec2f(0, 0)) const {

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
    inline Vec2f reflect (const Vec2f &n) const {
        double f = this->dot(n);
        return Vec2f(_x - n._x * 2.0 * f, _y - n._y * 2.0 * f);
    }

    // Length of the vector in units.
    //
    // @return length
    inline double length () const {
        return sqrt(_lengthSquared());
    }

    // Normalizes vector of the calling vector and returns it.
    // Does not change the calling vector.
    //
    // @return normalized vector, length 1
    inline Vec2f normalize () const {
        double l = length();
        if (l < std::numeric_limits<float>::epsilon()) {
            return *this;
        }
        return (*this / l);
    }

    // Scales the vector with a scalar s.
    //
    // @param scalar
    // @return scaled vector
    inline Vec2f operator* (float s) const {
        return Vec2f(_x * s, _y * s);
    }

    // Checks if this vector is laying on another.
    //
    // @param the other vector
    // @return  true if this vector lays on v
    //          false if this vector lays not on v
    inline bool on (Vec2f v) const {
        if (zero() && v.zero()) {
            return true;
        }
        else if (zero() || v.zero()) {
            return false;
        }
        double x_comp = _x / v.x(), y_comp = _y / v.y();
        Vec2f  x(_x, v.x()), y(_y, v.y());
        return (fabs(x_comp - y_comp) < 0.00001 || x.zero() || y.zero()) && (x_comp > 0 || x.zero()) &&
               (y_comp > 0 || y.zero()) && length() < v.length();
    }

    // Checks if this vector is laying in another.
    //
    // @param the other vector
    // @return  true if this vector lays in v
    //          false if this vector lays not in v
    inline bool in (const Vec2f& v) const {
        if (zero() && v.zero()) {
            return true;
        }
        else if (zero() || v.zero()) {
            return false;
        }

        double x_comp = _x / v.x(), y_comp = _y / v.y();
        Vec2f  x(_x, v.x()), y(_y, v.y());

        if ((y_comp < 0 && x_comp < 0) || (x.zero() && y_comp < 0) || (y.zero() && x_comp < 0)) {
            return false;
        }
        else if ((y_comp < 1 && x_comp < 1) || (x.zero() && y_comp < 1) || (y.zero() && x_comp < 1)) {
            return true;
        }
        else {
            return false;
        }

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
    inline double angle_signed (const Vec2f &v, bool sign_coord_x = true) const {
        int sign = -1;
        if (sign_coord_x && v._x >= _x) {
            sign = 1;
        }
        if (!sign_coord_x && v._y >= _y) {
            sign = 1;
        }

        return sign * angle(v);
    }

    // Calculates the angle between the calling vector and v.
    // Positive results mean clockwise,
    // negative results mean anticlockwise rotation.
    //
    // @param v: vector/point to find angle to.
    // @return angle in radians
    inline double real_angle (const Vec2f &v) const {
        return (atan2(_y, _x) - atan2(v.y(), v.x()));
    }

    // Calculates the angle between the calling vector and v.
    //
    // @param vector to find angle to.
    // @return angle in radians
    inline double angle (const Vec2f &v) const {
        double len_a = length();
        double len_b = v.length();

        if (len_a < 0.00001 || len_b < 0.00001) {
            return 0.0;
        }

        double angle = acos(dot(v) / (len_a * len_b));

        return angle < 0.00001 ? 0.0f : angle;
    }

    // Check equality of two vectors. Same functionality as operator==
    // but mainly used for testing.
    //
    // @param vector v to check equality against.
    // @return true if same values, else false
    inline bool equals (const Vec2f &v) const {
        return fabs(_x - v._x) < 0.00001 && fabs(_y - v._y) < 0.00001;
    }

    // Check equality of two vectors.
    //
    // @param vector v to check equality against.
    // @return true if same values, else false
    inline bool operator== (const Vec2f &b) const {
        return equals(b);
    }

    // Check inequality of two vectors.
    //
    // @param vector v to check equality against.
    // @return true if different values, else false
    inline bool operator!= (const Vec2f &b) const {
        return !equals(b);
    }

    // Implicit string conversion. Minimal string representation
    // of a vector. Same functionality as to_string().
    //
    // @return string representation of the vector.
    operator std::string () const;

    // Minimal string representation of a vector.
    //
    // @return string representation of the vector.
    std::string to_string () const;
};

// Vec3f represents a 3D vector in space which also can be a point.
// It implements basic arithmetic and helper functions, to create and manipulate
// points and vectors.
//
// For position information precision float was chosen because of its small memory footprint
// and therefore better caching in algorithms and use cases where millimeter precision
// is enough when 1 meter equals 1.0f.
//
// This class should never have any kind of setter methods. Algorithms and conversions need
// to be implemented inside the library.
class Vec3f {

protected:
    float _x;
    float _y;
    float _z;

private:
    // Squared length for utilization in multiple algorithms
    //
    // @return length squared
    inline double _lengthSquared () const {
        return this->dot(*this);
    }

public:
    // Constructor with initialization of both positions.
    //
    // @param x: position in x dimension
    // @param y: position in y dimension
    // @param z: position in z dimension
    Vec3f (float x, float y, float z);

    // Constructor initializes with x = y = z = 0.
    //
    Vec3f ();

    // Destructor - does nothing in this case.
    //
    ~Vec3f ();

    // Copy constructor initializes with x, y and z from v.
    //
    // @param vector to copy from
    Vec3f (const Vec3f &v);

    // Constructs a vector from a point to another point.
    //
    // @param from: origin
    // @param to: destination
    Vec3f (const Vec3f &from, const Vec3f &to);

    // Getter for x position.
    //
    // @return x position
    inline float x () const {
        return _x;
    }

    // Getter for y position.
    //
    // @return y position
    inline float y () const {
        return _y;
    }

    // Getter for z position.
    //
    // @return z position
    inline float z () const {
        return _z;
    }

    // Addition of x, y and z positions.
    //
    // @param other vector
    // @return new vector with the values from the addition
    inline Vec3f operator+ (const Vec3f &o) const {
        return Vec3f(_x + o._x, _y + o._y, _z + o._z);
    }

    // Addition with another vector but the calling class gets the result.
    //
    // @param other vector
    // @return *this for concatenation
    inline Vec3f &operator+= (const Vec3f &o) {
        _x += o._x;
        _y += o._y;
        _z += o._z;
        return *this;
    }

    // Assigns the calling class to the values of the given vector.
    //
    // @param other vector
    // @return *this for concatenation
    inline Vec3f &operator= (const Vec3f &o) {
        if (equals(o)) {
            return *this;
        }
        _x = o._x;
        _y = o._y;
        _z = o._z;
        return *this;
    }

    // Subtraction with another vector but the calling class gets the result.
    //
    // @param other vector
    // @return *this for concatenation
    inline Vec3f &operator-= (const Vec3f &o) {
        _x -= o._x;
        _y -= o._y;
        _z -= o._z;
        return *this;
    }

    // Subtraction with another vector.
    //
    // @param other vector
    // @return new vector with the result of the subtraction.
    inline Vec3f operator- (const Vec3f &o) const {
        return Vec3f(_x - o._x, _y - o._y, _z - o._z);
    }

    // Unary inversion. All x, y and z change the sign.
    //
    // @return new vector with inverted direction
    inline Vec3f operator- () const {
        return Vec3f(-_x, -_y, -_z);
    }

    // Subtraction with a scalar. All elements get subtracted by s.
    // Equal to subtracting by a Vector [s, s, s].
    //
    // @param scalar
    inline Vec3f operator- (float s) const {
        return Vec3f(_x - s, _y - s, _z - s);
    }

    // Division by a scalar. All elements get divided by s.
    //
    // @param scalar
    inline Vec3f operator/ (float s) const {
        return s < std::numeric_limits<float>::epsilon() ? Vec3f() : Vec3f(_x / s, _y / s, _z / s);
    }

    // Dot product with another vector. Mind the order. a.dot(b) == a . b.
    //
    // @param vector
    inline double dot (const Vec3f &v) const {
        return (_x * v._x) + (_y * v._y) + (_z * v._z);
    }

    // Checks whether a vector has values near zero (std::limits::epsilon).
    //
    // @return true if values are near zero, else false
    inline bool zero () const {
        return fabs(_x) < std::numeric_limits<float>::epsilon() &&
               fabs(_y) < std::numeric_limits<float>::epsilon() &&
               fabs(_z) < std::numeric_limits<float>::epsilon();
    }

    // Calculates the projected point from the calling vector on vector v.
    // This is also the closest point possible to the destination on the line that v describes.
    //
    // @param vector to project on
    // @return point on v
    inline Vec3f projected_point (const Vec3f &v) const {
        return v.normalize() * dot(v.normalize());
    }

    // Calculates the direct vector to any point on v.
    // The resulting vector has the smallest length possible.
    //
    // @param vector that describes the line
    // @return vector from destination of calling vector to v
    inline Vec3f closest_vec3_to (const Vec3f &v) const {
        return Vec3f(*this, projected_point(v));
    }

    // Rotates a vector around a pivot (origin is default), preserves height.
    // Positive radians mean clockwise,
    // negative radians mean anticlockwise rotation.
    //
    // @param radians to rotate
    // @param pivot to rotate around
    // @return rotated vector
    inline Vec3f rotate (float rad, const Vec3f &pivot = Vec3f(0, 0, 0)) const {

        double s = sin(rad);
        double c = cos(rad);

        // translate back to origin and rotate clockwise
        float x = c * (_x - pivot._x) + (_y - pivot._y) * s;
        float y = s * -(_x - pivot._x) + (_y - pivot._y) * c;

        // translate back
        return Vec3f(x + pivot._x, y + pivot._y, _z);
    }

    // Reflects a vector on another vector and returns the resulting bounce vector
    // in the normal direction of v.
    //
    // @param normal n
    // @return vector reflected on normal
    inline Vec3f reflect (const Vec3f &n) const {
        double f = this->dot(n);
        return Vec3f(_x - n._x * 2.0 * f, _y - n._y * 2.0 * f, _z - n._z * 2.0 * f);
    }

    // Length of the vector in units.
    //
    // @return length
    inline double length () const {
        return sqrt(_lengthSquared());
    }

    // Normalizes vector of the calling vector and returns it.
    // Does not change the calling vector.
    //
    // @return normalized vector, length 1
    inline Vec3f normalize () const {
        double l = length();
        if (l < std::numeric_limits<float>::epsilon()) {
            return *this;
        }
        return (*this / l);
    }

    // Scales the vector with a scalar s.
    //
    // @param scalar
    // @return scaled vector
    inline Vec3f operator* (float s) const {
        return Vec3f(_x * s, _y * s, _z * s);
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
    inline double angle_signed (const Vec3f &v, bool sign_coord_x = true) const {
        int sign = -1;
        if (sign_coord_x && v._x >= _x) {
            sign = 1;
        }
        if (!sign_coord_x && v._y >= _y) {
            sign = 1;
        }

        return sign * angle(v);
    }

    // Calculates the angle between the calling vector and v.
    //
    // @param vector to find angle to.
    // @return angle in radians
    inline double angle (const Vec3f &v) const {
        double len_a = length();
        double len_b = v.length();

        if (len_a < 0.00001 || len_b < 0.00001) {
            return 0.0;
        }

        double n     = dot(v) / (len_a * len_b);
        // branchless clamp
        double lower = -1, upper = 1;
        n = 0.5 * (n + lower + fabs(n - lower));
        n = 0.5 * (n + upper - fabs(upper - n));

        double angle = acos(n);

        return angle < 0.00001 ? 0.0f : angle;
    }

    // Check equality of two vectors. Same functionality as operator==
    // but mainly used for testing.
    //
    // @param vector v to check equality against.
    // @return true if same values, else false
    inline bool equals (const Vec3f &v) const {
        return fabs(_x - v._x) < 0.00001 && fabs(_y - v._y) < 0.00001 && fabs(_z - v._z) < 0.00001;
    }

    // Check equality of two vectors.
    //
    // @param vector v to check equality against.
    // @return true if same values, else false
    inline bool operator== (const Vec3f &b) const {
        return equals(b);
    }

    // Check inequality of two vectors.
    //
    // @param vector v to check equality against.
    // @return true if different values, else false
    inline bool operator!= (const Vec3f &b) const {
        return !equals(b);
    }

    // Implicit string conversion. Minimal string representation
    // of a vector. Same functionality as to_string().
    //
    // @return string representation of the vector.
    operator std::string () const;

    // Minimal string representation of a vector.
    //
    // @return string representation of the vector.
    std::string to_string () const;
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
class Polygon2 {
private:
    std::vector<Vec2f> vertices;
    vectortype         type = abs;

public:
    Polygon2 () = default;

    // Builds a polygon with the possibility to accept relative vectors.
    // O(1) time.
    //
    // @param polygontype what type of vectors this polygon should accept
    Polygon2 (vectortype polygontype);

    // Builds a polygon from all sides of a thing by combining the sides.
    // O(m) time.
    //
    // @param sides of a structure described in points
    Polygon2 (const std::vector<std::vector<Vec2f>> &sides);

    // Builds a polygon from all sides of a thing by combining the sides with the possibility to accept relative vectors.
    // O(m) time.
    //
    // @param sides of a structure described in points
    // @param polygontype what type of vectors this polygon should accept
    Polygon2 (const std::vector<std::vector<Vec2f>> &sides, vectortype polygontype);

    // Builds a polygon from points.
    // O(m) time.
    //
    // @param vector of all points.
    Polygon2 (const std::vector<Vec2f> &vertices);

    // Builds a polygon from points with the possibility to accept relative vectors.
    // O(m) time.
    //
    // @param vector of all points.
    // @param polygontype what type of vectors this polygon should accept
    Polygon2 (const std::vector<Vec2f> &vertices, vectortype polygontype);

    // Add a single vertex to the vertices programmatically.
    // O(1) time.
    //
    // @param vertex to add to the vertices.
    inline void add_vertex (const Vec2f &vertex) {
        vertices.push_back(vertices.size() ? vertices.back() * type + vertex : vertex);
    }

    // Checks whether a point is inside a polygon.
    // O(n) time.
    //
    // @param point to check
    // @return true if inside, else false
    bool contains (const Vec2f &p);

    // Smoothes the polygon with catmull-rom splines.
    // [WARNING] This function is not tested (because it is hard to do).
    // O(n*(1/distance)) time.
    //
    // @param alpha type of spline (0.0 = uniform, 0.5 = centripetal, 1.0 = chordal). Must be between 0.0 and 1.0
    // @param tension how tight the spline should be, 1.0 would result in linear splines. Must be between 0.0 and 1.0.
    // @param distance distance between new points in 1/100 of the distance between to original points
    void smooth (float alpha, float tension, unsigned int distance);

    // Adds a vertex to the polygon.
    // O(n) time complexity, because it returns a copy of the polygon.
    //
    // @param vertex to add
    // @return new polygon with added vertex
    inline Polygon2 operator+ (const Vec2f &vertex) const {
        Polygon2 ret = Polygon2(vertices);
        ret.add_vertex(vertices.size() ? vertices.back() * type + vertex : vertex);

        return ret;
    }

    // Adds a vertex to the polygon but inclusive and pretty.
    // O(1) time.
    //
    // @param vertex to add
    // @return new polygon with added vertex
    inline Polygon2 &operator<< (const Vec2f &vertex) {
        return (*this += vertex);
    }

    // Adds a vertex to the polygon but inclusive.
    // O(1) time.
    //
    // @param vertex to add
    // @return *this for concatenation
    inline Polygon2 &operator+= (const Vec2f &vertex) {
        add_vertex(vertices.size() ? vertices.back() * type + vertex : vertex);
        return *this;
    }

    // Moves all points of the Polygon.
    // O(n) time.
    //
    // @param movement of the polygon
    void move (const Vec2f &movement);

    // Rotates all points of the Polygon.
    // Positive radians mean clockwise,
    // negative radians mean anticlockwise rotation.
    // O(n) time.
    //
    // @param radians to rotate
    // @param pivot to rotate around
    void rotate (float rad, const Vec2f &pivot = Vec2f());

    // Sets the type of this polygon
    // O(1) time.
    //
    // @param polygontype the type the polygon should accept
    inline void setType (vectortype polygontype) {
        type = polygontype;
    }

    // Sizes the polygon along normals to allow for deadzones
    // Positive values mean enlargement.
    // Negative values mean shrinking.
    // O(n) time.
    //
    // @param dist how far the polygon should be scaled
    // @return the sized polygon
    [[deprecated("This function has been renamed to scaled()")]]
    Polygon2 sized (float dist);

    // Sizes the polygon along normals to allow for deadzones
    // Positive values mean enlargement.
    // Negative values mean shrinking.
    // O(n) time.
    //
    // @param dist how far the polygon should be scaled
    // @return the sized polygon
    Polygon2 scaled (float dist);

    // Minimal string representation of a polygon.
    // O(n) time.
    //
    // @return string representation of the polygon.
    std::string to_string () const;

    // Copying a polygon
    // O(n) time.
    //
    // @param clone the polygon to clone
    // @return this for concatenation
    Polygon2 &operator= (const Polygon2 &clone) {
        vertices = clone.vertices;
        type     = clone.type;
        return *this;
    }

    // Calculates the area of this polygon
    // O(n) time.
    //
    // @return the area covered by this polygon
    double area ();

    // Get the underlaying points
    // O(1) time.
    //
    // @return the vertices of the polygon
    inline std::vector<Vec2f> getVertices () const {
        return vertices;
    }

    // Copy constructor.
    // O(m) time.
    //
    // @param polygon to copy from.
    Polygon2 (const Polygon2 &v);
};


// Multiple points can be connected to define an object.
// The main structure describing that object is called a shape.
// [WARNING] Tis class is not yet tested
//
// Landau notation for all functions is given.
// Vector intern reallocation is ignored.
// n is scaling linearly with the amount of points in the polygon
// m is scaling linearly with the total amount of points passed
//
// !!! A shape only knows points, so each Vec3f will be interpreted as a point. !!!
// Points will be called vertex or vertices.
//
// It is assumed that each sub std::vector represents a layer of the shape at the same height
//
// The class manages the points and gives access to general manipulations of the structure.
// For specific functions, derive from this class.
class Shape3 {
private:
    std::vector<std::vector<Vec3f>> vertices;

public:
    Shape3 () = default;

    // Builds a polygon from all sides of a thing by combining the sides.
    // O(m) time.
    //
    // @param sides of a structure described in points
    Shape3 (const std::vector<std::vector<Vec3f>> &layers);

    // Add a single vertex to the vertices programmatically.
    // If no layer is given it is assumed to be added to the last layer
    // O(1) time.
    //
    // @param vertex to add to the vertices.
    inline void add_vertex (const Vec3f &vertex, int layer = -1) {
        layer = layer < 0 ? vertices.size() : layer;
        vertices[layer].push_back(vertex);
    }

    // Checks whether a point is inside a shape.
    // O(n) time.
    //
    // @param point to check
    // @return true if inside, else false
    // TODO interpolate between layers to find check as a polygon
    // bool contains(const Vec3f &p);

    // Smoothes the shape with catmull-rom splines.
    // O(n*(1/distance)) time.
    //
    // @param alpha type of spline (0.0 = uniform, 0.5 = centripetal, 1.0 = chordal). Must be between 0.0 and 1.0
    // @param tension how tight the spline should be, 1.0 would result in linear splines. Must be between 0.0 and 1.0.
    // @param distance distance between new points in 1/100 of the distance between to original points
    // TODO implement this on a layer by layer basis
    // void smooth(float alpha, float tension, unsigned int distance);

    // Adds a vertex to the shape on the last layer.
    // O(n) time complexity, because it returns a copy of the shape.
    //
    // @param vertex to add
    // @return new shape with added vertex
    inline Shape3 operator+ (const Vec3f &vertex) const {
        Shape3 ret = Shape3(vertices);
        ret.add_vertex(vertex);

        return ret;
    }

    // Adds a vertex to the shape on the last layer but inclusive and pretty.
    // O(1) time.
    //
    // @param vertex to add
    // @return new polygon with added vertex
    inline Shape3 &operator<< (const Vec3f &vertex) {
        return (*this += vertex);
    }

    // Adds a vertex to the shape on the last layer but inclusive.
    // O(1) time.
    //
    // @param vertex to add
    // @return *this for concatenation
    inline Shape3 &operator+= (const Vec3f &vertex) {
        add_vertex(vertex);
        return *this;
    }

    // Moves all points of the shape.
    // O(n) time.
    //
    // @param movement of the shape
    void move (const Vec3f &movement);

    // Rotates all points of the shape.
    // Positive radians mean clockwise,
    // negative radians mean anticlockwise rotation.
    // Height of the points is kept unaffected
    // O(n) time.
    //
    // @param radians to rotate
    // @param pivot to rotate around
    void rotate (float rad, const Vec3f &pivot = Vec3f());

    // Sizes the shape along normals to allow for deadzones
    // Positive values mean enlargement.
    // Negative values mean shrinking.
    // Height of the points is kept unaffected (for now)
    // O(n) time.
    //
    // @param dist how far the shape should be scaled
    // @return the sized shape
    Shape3 scaled (float dist);

    // Minimal string representation of a shape.
    // O(n) time.
    //
    // @return string representation of the shape.
    std::string to_string () const;

    // Copying a shape
    // O(n) time.
    //
    // @param clone the shape to clone
    // @return this for concatenation
    Shape3 &operator= (const Shape3 &clone) {
        vertices = clone.vertices;
        return *this;
    }

    // Calculates the volume of this shape
    // O(n) time.
    //
    // @return the volume covered by this shape
    // TODO implement this function
    // double volume();

    // Get the underlying points
    // O(1) time.
    //
    // @return the vertices of the shape
    inline std::vector<std::vector<Vec3f>> getVertices () const {
        return vertices;
    }

    // Copy constructor.
    // O(m) time.
    //
    // @param shape to copy from.
    Shape3 (const Shape3 &v);
};

} // namespace geo

#endif // geo_h