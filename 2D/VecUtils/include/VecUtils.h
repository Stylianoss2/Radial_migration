#ifndef VecUtils_h
#define VecUtils_h

#include <iostream>
#include <vector>
#include <cmath>

struct Vec3 
{
    double x, y, z;
};
Vec3 operator+(const Vec3 &v1, const Vec3 &v2) 
{
    return {v1.x + v2.x, v1.y + v2.y, 0.0};
}
Vec3 operator-(const Vec3 &v1, const Vec3 &v2) 
{
    return {v1.x - v2.x, v1.y - v2.y, 0.0};
}
Vec3 operator*(Vec3 const& v, double scalar)
{
    return {v.x * scalar, v.y * scalar, 0.0};
}
Vec3 operator/(Vec3 const& v, double scalar) 
{
    return {v.x / scalar, v.y / scalar, 0.0}; //Use 0.0 in 0 coordinate to avoid 0 / division by scalar
}

struct Cyl 
{
    double R, phi, z;
};
Cyl operator+(const Cyl &c1, const Cyl &c2)
{
    return {c1.R + c2.R, c1.phi + c2.phi, 0.0};
}
Cyl operator-(const Cyl &c1, const Cyl &c2) 
{
    return { c1.R - c2.R, c1.phi - c2.phi, 0.0};
}
Cyl operator*(const Cyl &c, double scalar)
{
    return {c.R * scalar, c.phi * scalar, 0.0};
}
Cyl operator/(const Cyl &c, double scalar) 
{
    return { c.R / scalar, c.phi / scalar, 0.0};
}

//////Conveersion between coordinate functions////////////
Vec3 cylindrical_to_cartesian(const Cyl &cyl) 
{
    double x = cyl.R * std::cos(cyl.phi);
    double y = cyl.R * std::sin(cyl.phi);
    double z = 0.0;

    return {x, y, z};
}

Cyl cartesian_to_cylindrical(const Vec3 &pos) 
{
    double R = std::sqrt(pos.x * pos.x + pos.y * pos.y);       
    double phi = std::atan2(pos.y, pos.x);                   
    if (phi < 0) phi += 2.0 * M_PI;  // φ = 0-2π
                                       

    return {R, phi, 0.0};
}

Vec3 cross_product(const Vec3& a, const Vec3& b)
{
    return 
    {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x 
    };
}
#endif