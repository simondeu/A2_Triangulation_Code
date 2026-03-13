#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>

using namespace easy3d;

// =====================================================================
// HELPER FUNCTION - must be ABOVE triangulation()
// =====================================================================
Matrix33 normalize(const std::vector<Vector2D>& pts, std::vector<Vector2D>& normalized_pts) {
    // compute centroid
    double cx = 0, cy = 0;
    for (auto& p : pts) {
        cx += p.x();
        cy += p.y();
    }
    cx /= pts.size();
    cy /= pts.size();

    // compute mean distance from centroid
    double mean_dist = 0;
    for (auto& p : pts)
        mean_dist += sqrt(pow(p.x()-cx, 2) + pow(p.y()-cy, 2));
    mean_dist /= pts.size();

    // scale factor
    double s = sqrt(2.0) / mean_dist;

    // build transform matrix
    Matrix33 T(s,  0,  -s*cx,
               0,  s,  -s*cy,
               0,  0,    1  );

    // apply to all points
    for (auto& p : pts) {
        Vector3D ph = T * p.homogeneous();
        normalized_pts.push_back(ph.cartesian());
    }

    return T;
}

// =====================================================================
// MAIN FUNCTION
// =====================================================================
bool Triangulation::triangulation(
        double fx, double fy,
        double cx, double cy,
        double s,
        const std::vector<Vector2D> &points_0,
        const std::vector<Vector2D> &points_1,
        std::vector<Vector3D> &points_3d,
        Matrix33 &R,
        Vector3D &t
) const
{
    // ------------------------------------------------------------------
    // 1. Input validation
    // ------------------------------------------------------------------
    if (points_0.size() < 8 || points_1.size() < 8) {
        std::cout << "Error: need at least 8 point pairs." << std::endl;
        return false;
    }
    if (points_0.size() != points_1.size()) {
        std::cout << "Error: point sets must be the same size." << std::endl;
        return false;
    }
    if (fx <= 0 || fy <= 0) {
        std::cout << "Error: focal lengths must be positive." << std::endl;
        return false;
    }

    // ------------------------------------------------------------------
    // 2. Normalize points
    // ------------------------------------------------------------------
    std::vector<Vector2D> norm_pts0, norm_pts1;
    Matrix33 T0 = normalize(points_0, norm_pts0);
    Matrix33 T1 = normalize(points_1, norm_pts1);

    std::cout << "=== Normalization ===" << std::endl;
    std::cout << "T0:\n" << T0 << std::endl;
    std::cout << "T1:\n" << T1 << std::endl;

    // ------------------------------------------------------------------
    // next steps go here (W matrix, SVD, F, E, R, t, triangulation)
    // ------------------------------------------------------------------
    Matrix WMatrix(points_0.size(), 9, 0.0);

    for (size_t i = 0; i < points_0.size(); ++i) {

    double u  = points_0[i][0];
    double v  = points_0[i][1];
    double up = points_1[i][0];
    double vp = points_1[i][1];


    WMatrix.set_row(i, {u * up,
        v * up,
        up,
        u * vp,
        v * vp,
        vp,
        u,
        v,
        1.0});
    };


    std::cout << "The W matrix" << WMatrix << std::endl;

    Matrix U(points_0.size(), points_0.size(), 0.0);
    Matrix S(points_0.size(), 9, 0.0);
    Matrix V(9, 9, 0.0);

    svd_decompose(WMatrix, U, S, V);


    return points_3d.size() > 0;
}