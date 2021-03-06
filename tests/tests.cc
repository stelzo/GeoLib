#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "geo.h"

#include <string>
#include <sstream>

using namespace geo;

TEST_CASE("rad degree conversion", "[base]")
{
    CHECK(deg2rad(1) == Approx(0.0174533));
    CHECK(deg2rad(400) == Approx(6.98132));

    CHECK(rad2deg(3) == Approx(171.887));
    CHECK(deg2rad(-1) == Approx(-0.0174533));
}

TEST_CASE("construct", "[vec2]")
{
    Vec2f a(2, 2);
    Vec2f c(3, 3);

    Vec2f g(1, 4);
    Vec2f f(2, 5);

    Vec2f b;
    Vec2f e(a);
    Vec2f d(a, c);

    CHECK(a.x() == Approx(2));
    CHECK(a.y() == Approx(2));

    CHECK(c.x() == Approx(3));
    CHECK(c.y() == Approx(3));

    CHECK(b.x() == Approx(0));
    CHECK(b.y() == Approx(0));

    CHECK(e.x() == Approx(2));
    CHECK(e.y() == Approx(2));

    CHECK(d.x() == Approx(1));
    CHECK(d.y() == Approx(1));

    CHECK(Vec2f(g, f).equals(Vec2f(1, 1)));

    CHECK_FALSE(a.zero());
    CHECK_FALSE(c.zero());

    CHECK(b.zero());
    CHECK_FALSE(e.zero());

    CHECK_FALSE(d.zero());

    auto y = Vec2f(1, 1);
    auto z = Vec2f(1, 1);

    CHECK(y == z);
    CHECK_FALSE(y == a);
    CHECK(y != a);

    std::string str_res = "[Vec2] x: 1 -- y: 4";

    CHECK(g.to_string() == str_res);

    std::stringstream ss;
    ss << (std::string)g;
    CHECK(ss.str() == str_res);
}

TEST_CASE("basic arithm", "[vec2]")
{

    Vec2f a(2, 2);
    Vec2f b(4, 5);

    // make sure the basics for the following tests are ok
    CHECK_FALSE(a.zero());
    CHECK_FALSE(b.zero());

    // tests depend on working equality operator
    Vec2f t(2, 2);
    CHECK(a.equals(t));

    SECTION("plus")
    {
        Vec2f res(6, 7);
        auto c = a + b;
        CHECK(c.equals(res));
    }

    SECTION("minus")
    {
        Vec2f res(2, 3);
        Vec2f c = b - a;
        CHECK(c.equals(res));
    }

    SECTION("scale (multiplication)")
    {
        Vec2f res(8, 10);
        Vec2f c = (b * 2);
        CHECK(c.equals(res));
    }

    SECTION("inline plus")
    {
        Vec2f res(6, 7);
        a += b;
        CHECK(a.equals(res));
    }

    SECTION("inline minus")
    {
        Vec2f res(2, 3);
        b -= a;
        CHECK(b.equals(res));
    }

    SECTION("unary minus -- invert")
    {
        Vec2f res(-4, -5);
        b = -b;
        CHECK(b.equals(res));
    }

    SECTION("minus scalar")
    {
        Vec2f res(2, 3);
        b = b - 2;
        CHECK(b.equals(res));
    }

    SECTION("divide scalar")
    {
        Vec2f res(3, 4);
        a = Vec2f(9, 12) / 3;
        CHECK(a.equals(res));
    }
}

TEST_CASE("functions", "[vec2]")
{

    Vec2f a(2, 2);
    Vec2f b(4, 4);

    SECTION("dot product easy")
    {
        CHECK(a.dot(b) == Approx(16));
    }

    SECTION("dot product harder")
    {
        CHECK(Vec2f(3, 5).dot(Vec2f(8, 2)) == Approx(34));
    }

    SECTION("rotate clockwise around point")
    {
        auto c = a.rotate(M_PI / 2, b);
        Vec2f res(2, 6);
        CHECK(c.equals(res));
    }

    SECTION("rotate counter-clockwise around origin")
    {
        auto c = a.rotate(-M_PI / 2);
        Vec2f res(-2, 2);
        CHECK(c.equals(res));
    }

    SECTION("length")
    {
        CHECK(a.length() == Approx(2.82843));
    }

    SECTION("normalize")
    {
        CHECK(a.normalize().length() == Approx(1));
    }

    SECTION("angle zero")
    {
        CHECK(a.angle(b) == Approx(0));
    }

    SECTION("angle vector")
    {
        CHECK(a.angle(Vec2f(2, 4)) == Approx(0.321751));
    }

    SECTION("real angle vector")
    {
        CHECK(Vec2f(1, 0).real_angle(Vec2f(0, 1)) == Approx(-M_PI_2));
        CHECK(Vec2f(1, 0).real_angle(Vec2f(0, -1)) == Approx(M_PI_2));
        CHECK(a.real_angle(Vec2f(2, 4)) == Approx(-0.321751));
        CHECK(a.real_angle(Vec2f(4, 2)) == Approx(0.321751));
    }

    SECTION("on")
    {
        CHECK(Vec2f(1, 0).on(Vec2f(2, 0)));
        CHECK(Vec2f(1, 1).on(Vec2f(2, 2)));
        CHECK_FALSE(Vec2f(1, 0).on(Vec2f(2, 2)));
        CHECK_FALSE(Vec2f(-1, -1).on(Vec2f(2, 2)));
        CHECK(a.on(Vec2f(4, 4)));
        CHECK_FALSE(a.on(Vec2f(4, 2)));
    }

    SECTION("create on")
    {
        Vec2f f = Vec2f(1, 1).closest_connection(Vec2f(2, 0));
        Vec2f g = Vec2f(3, 1).closest_connection(Vec2f(2, 0));
        Vec2f h = Vec2f(-1, 1).closest_connection(Vec2f(2, 0));
        Vec2f res_f(0, -1);
        Vec2f res_g(-1, -1);
        Vec2f res_h(1, -1);
        CHECK(f.equals(res_f));
        CHECK(g.equals(res_g));
        CHECK(h.equals(res_h));
    }

    SECTION("reflection")
    {
        auto c = -Vec2f(-1, 1).reflect(Vec2f(0, 1));
        CHECK(c.equals(Vec2f(1, 1)));
    }

    SECTION("projected point")
    {
        CHECK(a.projected_point(Vec2f(2, 4)).equals(Vec2f(1.2, 2.4)));
    }
}

TEST_CASE("construct3", "[vec3]")
{
    Vec3f a(2, 2, 2);
    Vec3f c(3, 3, 3);

    Vec3f g(1, 4, 7);
    Vec3f f(2, 5, 8);

    Vec3f b;
    Vec3f e(a);
    Vec3f d(a, c);

    CHECK(a.x() == Approx(2));
    CHECK(a.y() == Approx(2));
    CHECK(a.z() == Approx(2));

    CHECK(c.x() == Approx(3));
    CHECK(c.y() == Approx(3));
    CHECK(c.z() == Approx(3));

    CHECK(b.x() == Approx(0));
    CHECK(b.y() == Approx(0));
    CHECK(b.z() == Approx(0));

    CHECK(e.x() == Approx(2));
    CHECK(e.y() == Approx(2));
    CHECK(e.z() == Approx(2));

    CHECK(d.x() == Approx(1));
    CHECK(d.y() == Approx(1));
    CHECK(d.z() == Approx(1));

    CHECK(Vec3f(g, f).equals(Vec3f(1, 1, 1)));

    CHECK_FALSE(a.zero());
    CHECK_FALSE(c.zero());

    CHECK(b.zero());
    CHECK_FALSE(e.zero());

    CHECK_FALSE(d.zero());

    auto y = Vec3f(1, 1, 1);
    auto z = Vec3f(1, 1, 1);

    CHECK(y == z);
    CHECK_FALSE(y == a);
    CHECK(y != a);

    std::string str_res = "[Vec3] x: 1 -- y: 4 -- z: 7";

    CHECK(g.to_string() == str_res);

    std::stringstream ss;
    ss << (std::string)g;
    CHECK(ss.str() == str_res);
}

TEST_CASE("basic arithm3", "[vec3]")
{

    Vec3f a(2, 2, 2);
    Vec3f b(4, 5, 6);

    // make sure the basics for the following tests are ok
    CHECK_FALSE(a.zero());
    CHECK_FALSE(b.zero());

    // tests depend on working equality operator
    Vec3f t(2, 2, 2);
    CHECK(a.equals(t));

    SECTION("plus")
    {
        Vec3f res(6, 7, 8);
        auto c = a + b;
        CHECK(c.equals(res));
    }

    SECTION("minus")
    {
        Vec3f res(2, 3, 4);
        Vec3f c = b - a;
        CHECK(c.equals(res));
    }

    SECTION("scale (multiplication)")
    {
        Vec3f res(8, 10, 12);
        Vec3f c = (b * 2);
        CHECK(c.equals(res));
    }

    SECTION("inline plus")
    {
        Vec3f res(6, 7, 8);
        a += b;
        CHECK(a.equals(res));
    }

    SECTION("inline minus")
    {
        Vec3f res(2, 3, 4);
        b -= a;
        CHECK(b.equals(res));
    }

    SECTION("unary minus -- invert")
    {
        Vec3f res(-4, -5, -6);
        b = -b;
        CHECK(b.equals(res));
    }

    SECTION("minus scalar")
    {
        Vec3f res(2, 3, 4);
        b = b - 2;
        CHECK(b.equals(res));
    }

    SECTION("divide scalar")
    {
        Vec3f res(3, 4, 5);
        a = Vec3f(9, 12, 15) / 3;
        CHECK(a.equals(res));
    }
}

TEST_CASE("functions3", "[vec3]")
{

    Vec3f a(2, 2, 2);
    Vec3f b(4, 4, 4);

    SECTION("dot product easy")
    {
        CHECK(a.dot(b) == Approx(24));
    }

    SECTION("dot product harder")
    {
        CHECK(Vec3f(3, 5, 7).dot(Vec3f(8, 2, -4)) == Approx(6));
    }

    SECTION("rotate clockwise around point")
    {
        auto c = a.rotate(M_PI / 2, b);
        Vec3f res(2, 6, 2);
        CHECK(c.equals(res));
    }

    SECTION("rotate counter-clockwise around origin")
    {
        auto c = a.rotate(-M_PI / 2);
        Vec3f res(-2, 2, 2);
        CHECK(c.equals(res));
    }

    SECTION("length")
    {
        CHECK(a.length() == Approx(3.464101615));
    }

    SECTION("normalize")
    {
        CHECK(a.normalize().length() == Approx(1));
    }

    SECTION("angle zero")
    {
        CHECK(a.angle(b) == Approx(0));
    }

    SECTION("angle vector")
    {
        CHECK(a.angle(Vec3f(2, 4, 2)) == Approx(0.3398369095));
    }

    SECTION("reflection")
    {
        auto c = -Vec3f(-1, 1, 0).reflect(Vec3f(0, 1, 0));
        CHECK(c.equals(Vec3f(1, 1, 0)));
    }

    SECTION("projected point")
    {
        CHECK(a.projected_point(Vec3f(2, 4, 0)).equals(Vec3f(1.2, 2.4, 0)));
    }
}

TEST_CASE("all", "[Polygon2]")
{
    std::vector<Vec2f> ps;
    ps.emplace_back(1, 2);
    ps.emplace_back(3, 4);
    ps.emplace_back(2, 1);

    // operators to add vertex
    Polygon2 cc;
    cc += Vec2f(1, 2);
    cc += Vec2f(3, 4);

    Polygon2 dd = cc + Vec2f(2, 1);
    Polygon2 ee = cc << Vec2f(2, 1);

    cc += Vec2f(2, 1);

    std::vector<std::vector<Vec2f>> ps_vec;
    ps_vec.push_back(ps);

    Polygon2 pol(ps);
    Polygon2 pol_v(ps_vec);

    SECTION("close outside")
    {
        CHECK_FALSE(pol.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(pol_v.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(cc.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(dd.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(ee.contains(Vec2f(2.5, 2.5)));
    }

    SECTION("easy inside")
    {
        CHECK(pol.contains(Vec2f(2, 2.5)));
        CHECK(pol_v.contains(Vec2f(2, 2.5)));
        CHECK(cc.contains(Vec2f(2, 2.5)));
        CHECK(dd.contains(Vec2f(2, 2.5)));
        CHECK(ee.contains(Vec2f(2, 2.5)));
    }

    Polygon2 square;
    square << Vec2f(0, 0);
    square << Vec2f(1, 0);
    square << Vec2f(1, 1);
    square << Vec2f(0, 1);

    Polygon2 mov(square), rot(square), rotoffcenter(square);
    mov.move(Vec2f(2, 2));
    rot.rotate(M_PI);
    rotoffcenter.rotate(M_PI, Vec2f(2.5, 2.5));

    SECTION("placement")
    {
        CHECK_FALSE(square.contains(Vec2f(-0.5, -0.5)));
        CHECK(square.contains(Vec2f(0.5, 0.5)));
        CHECK_FALSE(square.contains(Vec2f(1.5, 1.5)));
        CHECK_FALSE(square.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(square.contains(Vec2f(4.5, 4.5)));
    }

    SECTION("movement")
    {
        CHECK_FALSE(mov.contains(Vec2f(-0.5, -0.5)));
        CHECK_FALSE(mov.contains(Vec2f(0.5, 0.5)));
        CHECK_FALSE(mov.contains(Vec2f(1.5, 1.5)));
        CHECK(mov.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(mov.contains(Vec2f(4.5, 4.5)));
    }

    SECTION("rotation")
    {
        CHECK(rot.contains(Vec2f(-0.5, -0.5)));
        CHECK_FALSE(rot.contains(Vec2f(0.5, 0.5)));
        CHECK_FALSE(rot.contains(Vec2f(1.5, 1.5)));
        CHECK_FALSE(rot.contains(Vec2f(2.5, 2.5)));
        CHECK_FALSE(rot.contains(Vec2f(4.5, 4.5)));
    }

    SECTION("rotation off center")
    {
        CHECK_FALSE(rotoffcenter.contains(Vec2f(-0.5, -0.5)));
        CHECK_FALSE(rotoffcenter.contains(Vec2f(0.5, 0.5)));
        CHECK_FALSE(rotoffcenter.contains(Vec2f(1.5, 1.5)));
        CHECK_FALSE(rotoffcenter.contains(Vec2f(2.5, 2.5)));
        CHECK(rotoffcenter.contains(Vec2f(4.5, 4.5)));
    }

    Polygon2 scaled = square.scaled(1.0);

    SECTION("sizing")
    {
        CHECK(scaled.contains(Vec2f(-0.5, -0.5)));
        CHECK(scaled.contains(Vec2f(-0.5, 1.5)));
        CHECK(scaled.contains(Vec2f(1.5, -0.5)));
        CHECK(scaled.contains(Vec2f(1.5, 1.5)));

        CHECK_FALSE(square.contains(Vec2f(-0.5, -0.5)));
        CHECK_FALSE(square.contains(Vec2f(-0.5, 1.5)));
        CHECK_FALSE(square.contains(Vec2f(1.5, -0.5)));
        CHECK_FALSE(square.contains(Vec2f(1.5, 1.5)));
    }

    SECTION("areas")
    {
        CHECK(square.area() == Approx(1));
        CHECK(mov.area() == Approx(1));
        CHECK(rot.area() == Approx(1));
        CHECK(rotoffcenter.area() == Approx(1));
    }
}