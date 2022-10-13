#pragma once

#include <math.h>
#include <assert.h>
#include <limits>
#include <algorithm>

// Helper functions ...

double pseudoAngle(double dx, double dy)
{
    double p = dx / (abs(dx) + abs(dy));
    return (dy > 0 ? 3 - p : 1 + p) / 4; // [0..1]
}

inline double dist(double ax, double ay, double bx, double by)
{
    double dx = ax - bx;
    double dy = ay - by;
    return dx * dx + dy * dy;
}

inline double orient2d(double ax, double ay, double bx, double by, double cx, double cy)
{
    return (ay - cy) * (bx - cx) - (ax - cx) * (by - cy);
}

bool inCircle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py)
{
    double dx = ax - px;
    double dy = ay - py;
    double ex = bx - px;
    double ey = by - py;
    double fx = cx - px;
    double fy = cy - py;

    double ap = dx * dx + dy * dy;
    double bp = ex * ex + ey * ey;
    double cp = fx * fx + fy * fy;

    return dx * (ey * cp - bp * fy) -
        dy * (ex * cp - bp * fx) +
        ap * (ex * fy - ey * fx) < 0;
}

double circumradius(double ax, double ay, double bx, double by, double cx, double cy)
{
    double dx = bx - ax;
    double dy = by - ay;
    double ex = cx - ax;
    double ey = cy - ay;

    double bl = dx * dx + dy * dy;
    double cl = ex * ex + ey * ey;
    double d = 0.5 / (dx * ey - dy * ex);

    double x = (ey * bl - dy * cl) * d;
    double y = (dx * cl - ex * bl) * d;

    return x * x + y * y;
}

void circumcenter(double ax, double ay, double bx, double by, double cx, double cy, double* x, double* y)
{
    double dx = bx - ax;
    double dy = by - ay;
    double ex = cx - ax;
    double ey = cy - ay;

    double bl = dx * dx + dy * dy;
    double cl = ex * ex + ey * ey;
    double d = 0.5 / (dx * ey - dy * ex);

    *x = ax + (ey * bl - dy * cl) * d;
    *y = ay + (dx * cl - ex * bl) * d;
}

inline void swap(int* arr, int i, int j)
{
    const int tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

void quicksort(int* ids, double* dists, int left, int right)
{
    if (right - left <= 10)
    {
        for (int i = left + 1; i <= right; i++)
        {
            int temp = ids[i];
            double tempDist = dists[temp];
            int  j = i - 1;
            while (j >= left && dists[ids[j]] > tempDist)
            {
                ids[j + 1] = ids[j];
                j--;
            }
            ids[j + 1] = temp;
        }
    }
    else
    {
        int median = (left + right) >> 1;
        int i = left + 1;
        int j = right;
        swap(ids, median, i);
        if (dists[ids[left]] > dists[ids[right]]) swap(ids, left, right);
        if (dists[ids[i]] > dists[ids[right]]) swap(ids, i, right);
        if (dists[ids[left]] > dists[ids[i]]) swap(ids, left, i);

        int temp = ids[i];
        double tempDist = dists[temp];
        while (true) {
            do i++; while (dists[ids[i]] < tempDist);
            do j--; while (dists[ids[j]] > tempDist);
            if (j < i) break;
            swap(ids, i, j);
        }
        ids[left + 1] = ids[j];
        ids[j] = temp;

        if (right - i + 1 >= j - left)
        {
            quicksort(ids, dists, i, right);
            quicksort(ids, dists, left, j - 1);
        }
        else
        {
            quicksort(ids, dists, left, j - 1);
            quicksort(ids, dists, i, right);
        }
    }
}

class Delaunator
{
private:
    const double EPSILON = pow(2, -52);
    const double Infinity = std::numeric_limits<double>::infinity();

    int EDGE_STACK[512];
    double* coords;
    int n;
    int hashSize;
    int hullStart;
    double cx, cy;

    int* hullPrev;
    int* hullNext;
    int* hullTri;
    int* hullHash;
    int* ids;

public:
    int trianglesLen;
    int* triangles;
    int* halfedges;

    ~Delaunator()
    {
        delete[] triangles;
        delete[] halfedges;
        delete[] hullPrev;
        delete[] hullNext;
        delete[] hullTri;
        delete[] hullHash;
        delete[] ids;
    }

    Delaunator(double* coords, int n)
    {
        this->coords = coords;
        this->n = n;

        hashSize = (int)ceil(sqrt(n));
        const int maxTriangles = std::max(2 * n - 5, 0);
        
        triangles = new int[maxTriangles * 3];
        halfedges = new int[maxTriangles * 3];
        hullPrev = new int[n];
        hullNext = new int[n];
        hullTri = new int[n]; 
        hullHash = new int[hashSize]; 
        ids = new int[n];

        // populate an array of point indices; calculate input data bbox
        double minX = Infinity;
        double minY = Infinity;
        double maxX = -Infinity;
        double maxY = -Infinity;

        for (int i = 0; i < n; i++)
        {
            double x = coords[2 * i];
            double y = coords[2 * i + 1];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
            ids[i] = i;
        }
        double cx0 = (minX + maxX) / 2;
        double cy0 = (minY + maxY) / 2;

        double minDist = Infinity;
        int i0, i1, i2;

        // pick a seed point close to the center
        for (int i = 0; i < n; i++) {
            double d = dist(cx0, cy0, coords[2 * i], coords[2 * i + 1]);
            if (d < minDist) {
                i0 = i;
                minDist = d;
            }
        }
        double i0x = coords[2 * i0];
        double i0y = coords[2 * i0 + 1];

        minDist = Infinity;

        // find the point closest to the seed
        for (int i = 0; i < n; i++)
        {
            if (i == i0) continue;
            double d = dist(i0x, i0y, coords[2 * i], coords[2 * i + 1]);
            if (d < minDist && d > 0)
            {
                i1 = i;
                minDist = d;
            }
        }
        double i1x = coords[2 * i1];
        double i1y = coords[2 * i1 + 1];

        double minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (int i = 0; i < n; i++)
        {
            if (i == i0 || i == i1) continue;
            double r = circumradius(i0x, i0y, i1x, i1y, coords[2 * i], coords[2 * i + 1]);
            if (r < minRadius)
            {
                i2 = i;
                minRadius = r;
            }
        }
        double i2x = coords[2 * i2];
        double i2y = coords[2 * i2 + 1];

        assert(minRadius != Infinity);

        // swap the order of the seed points for counter-clockwise orientation
        if (orient2d(i0x, i0y, i1x, i1y, i2x, i2y) < 0)
        {
            int i = i1;
            int x = (int)i1x;
            int y = (int)i1y;
            i1 = i2;
            i1x = i2x;
            i1y = i2y;
            i2 = i;
            i2x = x;
            i2y = y;
        }

        circumcenter(i0x, i0y, i1x, i1y, i2x, i2y, &cx, &cy);

        double* dists = new double[n];

        for (int i = 0; i < n; i++)
            dists[i] = dist(coords[2 * i], coords[2 * i + 1], cx, cy);

        // sort the points by distance from the seed triangle circumcenter
        //std::sort(ids, ids+n, [dists](const int& a, const int& b) -> bool { return dists[a] <= dists[b]; });
        quicksort(ids, dists, 0, n - 1);

        delete[] dists;

        // set up the seed triangle as the starting hull
        hullStart = i0;
        int hullSize = 3;

        hullNext[i0] = hullPrev[i2] = i1;
        hullNext[i1] = hullPrev[i0] = i2;
        hullNext[i2] = hullPrev[i1] = i0;

        hullTri[i0] = 0;
        hullTri[i1] = 1;
        hullTri[i2] = 2;

        //hullHash.fill(-1);
        for (int i = 0; i < hashSize; i++)
            hullHash[i] = -1;

        hullHash[hashKey(i0x, i0y)] = i0;
        hullHash[hashKey(i1x, i1y)] = i1;
        hullHash[hashKey(i2x, i2y)] = i2;

        trianglesLen = 0;
        addTriangle(i0, i1, i2, -1, -1, -1);

        double xp, yp;
        for (int k = 0; k < n; k++)
        {
            int i = ids[k];
            double x = coords[2 * i];
            double y = coords[2 * i + 1];

            // skip near-duplicate points
            if (k > 0 && abs(x - xp) <= EPSILON && abs(y - yp) <= EPSILON) 
                continue;
            xp = x;
            yp = y;

            // skip seed triangle points
            if (i == i0 || i == i1 || i == i2) continue;

            // find a visible edge on the convex hull using edge hash
            int start = 0;
            for (int j = 0, key = hashKey(x, y); j < hashSize; j++) {
                start = hullHash[(key + j) % hashSize];
                if (start != -1 && start != hullNext[start]) break;
            }

            start = hullPrev[start];
            int e = start;
            int q;
            while (q = hullNext[e], orient2d(x, y, coords[2 * e], coords[2 * e + 1], coords[2 * q], coords[2 * q + 1]) >= 0)
            {
                e = q;
                if (e == start)
                {
                    e = -1;
                    break;
                }
            }
            if (e == -1) continue; // likely a near-duplicate point; skip it

            // add the first triangle from the point
            int t = addTriangle(e, i, hullNext[e], -1, -1, hullTri[e]);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            hullTri[i] = legalize(t + 2);
            hullTri[e] = t; // keep track of boundary triangles on the hull
            hullSize++;

            // walk forward through the hull, adding more triangles and flipping recursively
            int n = hullNext[e];
            while (q = hullNext[n], orient2d(x, y, coords[2 * n], coords[2 * n + 1], coords[2 * q], coords[2 * q + 1]) < 0)
            {
                t = addTriangle(n, i, q, hullTri[i], -1, hullTri[n]);
                hullTri[i] = legalize(t + 2);
                hullNext[n] = n; // mark as removed
                hullSize--;
                n = q;
            }

            // walk backward from the other side, adding more triangles and flipping
            if (e == start)
            {
                while (q = hullPrev[e], orient2d(x, y, coords[2 * q], coords[2 * q + 1], coords[2 * e], coords[2 * e + 1]) < 0)
                {
                    t = addTriangle(q, i, e, -1, hullTri[e], hullTri[q]);
                    legalize(t + 2);
                    hullTri[q] = t;
                    hullNext[e] = e; // mark as removed
                    hullSize--;
                    e = q;
                }
            }

            // update the hull indices
            hullStart = hullPrev[i] = e;
            hullNext[e] = hullPrev[n] = i;
            hullNext[i] = n;

            // save the two new edges in the hash table
            hullHash[hashKey(x, y)] = i;
            hullHash[hashKey(coords[2 * e], coords[2 * e + 1])] = e;
        }
    }

    int hashKey(double x, double y)
    {
        return (int)(pseudoAngle(x - cx, y - cy) * hashSize) % hashSize;
    }

    int legalize(int a)
    {
        int i = 0;
        int ar = 0;

        // recursion eliminated with a fixed-size stack
        while (true)
        {
            int b = halfedges[a];
            int a0 = a - a % 3;
            ar = a0 + (a + 2) % 3;

            if (b == -1)
            {
                if (i == 0) break;
                a = EDGE_STACK[--i];
                continue;
            }

            int b0 = b - b % 3;
            int al = a0 + (a + 1) % 3;
            int bl = b0 + (b + 2) % 3;

            int p0 = triangles[ar];
            int pr = triangles[a];
            int pl = triangles[al];
            int p1 = triangles[bl];

            bool illegal = inCircle(
                coords[2 * p0], coords[2 * p0 + 1],
                coords[2 * pr], coords[2 * pr + 1],
                coords[2 * pl], coords[2 * pl + 1],
                coords[2 * p1], coords[2 * p1 + 1]);

            if (illegal)
            {
                triangles[a] = p1;
                triangles[b] = p0;

                int hbl = halfedges[bl];

                // edge swapped on the other side of the hull (rare); fix the halfedge reference
                if (hbl == -1)
                {
                    int e = hullStart;
                    do {
                        if (hullTri[e] == bl) {
                            hullTri[e] = a;
                            break;
                        }
                        e = hullPrev[e];
                    } while (e != hullStart);
                }
                link(a, hbl);
                link(b, halfedges[ar]);
                link(ar, bl);

                int br = b0 + (b + 1) % 3;

                // don't worry about hitting the cap: it can only happen on extremely degenerate input
                assert(i < 512);
                EDGE_STACK[i++] = br;
            }
            else
            {
                if (i == 0) break;
                a = EDGE_STACK[--i];
            }
        }

        return ar;
    }

    void link(int a, int b)
    {
        halfedges[a] = b;
        if (b != -1) halfedges[b] = a;
    }

    // add a new triangle given vertex indices and adjacent half-edge ids
    int addTriangle(int i0, int i1, int i2, int a, int b, int c)
    {
        int t = trianglesLen;

        triangles[t] = i0;
        triangles[t + 1] = i1;
        triangles[t + 2] = i2;

        link(t, a);
        link(t + 1, b);
        link(t + 2, c);

        trianglesLen += 3;

        return t;
    }
};

