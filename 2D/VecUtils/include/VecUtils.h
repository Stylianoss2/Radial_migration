#ifndef VecUtils_h
#define VecUtils_h

#include <iostream>
#include <vector>
#include <cmath>

struct Vec2 
{
    double x, y;
};
inline Vec2 operator+(const Vec2 &v1, const Vec2 &v2) 
{
    return {v1.x + v2.x, v1.y + v2.y};
}
inline Vec2 operator-(const Vec2 &v1, const Vec2 &v2) 
{
    return {v1.x - v2.x, v1.y - v2.y};
}
inline Vec2 operator*(Vec2 const& v, double scalar)
{
    return {v.x * scalar, v.y * scalar};
}
inline Vec2 operator/(Vec2 const& v, double scalar) 
{
    return {v.x / scalar, v.y / scalar}; 
}

struct Cyl 
{
    double R, phi;
};
inline Cyl operator+(const Cyl &c1, const Cyl &c2)
{
    return {c1.R + c2.R, c1.phi + c2.phi};
}
inline Cyl operator-(const Cyl &c1, const Cyl &c2) 
{
    return {c1.R - c2.R, c1.phi - c2.phi};
}
inline Cyl operator*(const Cyl &c, double scalar)
{
    return {c.R * scalar, c.phi * scalar};
}
inline Cyl operator/(const Cyl &c, double scalar) 
{
    return { c.R / scalar, c.phi / scalar};
}

//////Conveersion between coordinate functions////////////
inline Vec2 cylindrical_to_cartesian(const Cyl &cyl) 
{
    double x = cyl.R * std::cos(cyl.phi);
    double y = cyl.R * std::sin(cyl.phi);
    return {x, y};
}

inline Cyl cartesian_to_cylindrical(const Vec2 &pos) 
{
    double R = std::hypot(pos.x, pos.y);  
    double phi = std::atan2(pos.y, pos.x);

    if (phi < 0.0) 
        phi += 2.0 * M_PI;

    return {R, phi};
}

inline Cyl cartesian_to_cylindrical_velocity(const Vec2& pos, const Vec2& vel)
{
    double R = std::sqrt(pos.x * pos.x + pos.y * pos.y);
    double phi = std::atan2(pos.y, pos.x);

    double v_R = (pos.x * vel.x + pos.y * vel.y) / R;
    double v_phi = (pos.x * vel.y - pos.y * vel.x) / R; 

    return {v_R, v_phi}; 
}
#endif