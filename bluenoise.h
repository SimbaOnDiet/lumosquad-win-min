#pragma once
#include "stdafx.h"
#include "RNG.h"
#include <cmath>
#include <vector>

#define kMaxPointsPerCell 9

class RangeList;
class ScallopedRegion;

class Vec2 {
public:
	Vec2() {};
	Vec2(float _x, float _y) : x(_x), y(_y) {};

	float x, y;

	float length() { return sqrt(x*x + y*y); }

	bool operator ==(const Vec2 &b) const { return x == b.x && y == b.y; }
	Vec2 operator +(Vec2 b) { return Vec2(x + b.x, y + b.y); }
	Vec2 operator -(Vec2 b) { return Vec2(x - b.x, y - b.y); }
	Vec2 operator *(Vec2 b) { return Vec2(x*b.x, y*b.y); }
	Vec2 operator /(Vec2 b) { return Vec2(x / b.x, y*b.y); }

	Vec2 operator +(float n) { return Vec2(x + n, y + n); }
	Vec2 operator -(float n) { return Vec2(x - n, y - n); }
	Vec2 operator *(float n) { return Vec2(x*n, y*n); }
	Vec2 operator /(float n) { return Vec2(x / n, y*n); }

	Vec2 &operator +=(Vec2 b) { x += b.x; y += b.y; return *this; }
	Vec2 &operator -=(Vec2 b) { x -= b.x; y -= b.y; return *this; }
	Vec2 &operator *=(Vec2 b) { x *= b.x; y *= b.y; return *this; }
	Vec2 &operator /=(Vec2 b) { x /= b.x; y /= b.y; return *this; }

	Vec2 &operator +=(float n) { x += n; y += n; return *this; }
	Vec2 &operator -=(float n) { x -= n; y -= n; return *this; }
	Vec2 &operator *=(float n) { x *= n; y *= n; return *this; }
	Vec2 &operator /=(float n) { x /= n; y /= n; return *this; }
};

/// \brief Daniel Dunbar's blue noise generator
/// 
/// The original code has been modified so that the 'boundary sampling'
/// method is the only one available.
class BLUE_NOISE {
protected:
	RNG m_rng;
	std::vector<int> m_neighbors;

	int(*m_grid)[kMaxPointsPerCell];
	int m_gridSize;
	float m_gridCellSize;

public:
	std::vector<Vec2> points;
	float radius;
	bool isTiled;

public:
	BLUE_NOISE(float radius, bool isTiled = true, bool usesGrid = true);
	virtual ~BLUE_NOISE() { };

	//

	bool pointInDomain(Vec2 &a);

	// return shortest distance between _a_ 
	// and _b_ (accounting for tiling)
	float getDistanceSquared(Vec2 &a, Vec2 &b) { Vec2 v = getTiled(b - a); return v.x*v.x + v.y*v.y; }
	float getDistance(Vec2 &a, Vec2 &b) { return sqrt(getDistanceSquared(a, b)); }

	// generate a random point in square
	Vec2 randomPoint();

	// return tiled coordinates of _v_
	Vec2 getTiled(Vec2 v);

	// return grid x,y for point
	void getGridXY(Vec2 &v, int *gx_out, int *gy_out);

	// add _pt_ to point list and grid
	void addPoint(Vec2 pt);

	// populate m_neighbors with list of
	// all points within _radius_ of _pt_
	// and return number of such points
	int findNeighbors(Vec2 &pt, float radius);

	// return distance to closest neighbor within _radius_
	float findClosestNeighbor(Vec2 &pt, float radius);

	// find available angle ranges on boundary for candidate 
	// by subtracting occluded neighbor ranges from _rl_
	void findNeighborRanges(int index, RangeList &rl);

	// extend point set by boundary sampling until domain is 
	// full
	void maximize();

	// apply one step of Lloyd relaxation
	void relax();

	void complete();

	void writeToBool(bool* noise, int size);
};
