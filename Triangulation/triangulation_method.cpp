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
    Matrix WMatrix(norm_pts0.size(), 9, 0.0);

    for (size_t i = 0; i < norm_pts0.size(); ++i) {

    double u  = norm_pts0[i][0];
    double v  = norm_pts0[i][1];
    double up = norm_pts1[i][0];
    double vp = norm_pts1[i][1];


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


    std::cout << "The W matrix: " << WMatrix << std::endl;

    Matrix U(norm_pts0.size(), points_0.size(), 0.0);
    Matrix S(norm_pts0.size(), 9, 0.0);
    Matrix V(9, 9, 0.0);

    svd_decompose(WMatrix, U, S, V);

<<<<<<< HEAD
    
=======
    Matrix F(3, 3, V.get_column(8).data());
    
    std::cout << "The F matrix: " << F << std::endl;

>>>>>>> f49071ac902ae3b36382e6300e247beb060b1ece
    return points_3d.size() > 0;
}

// TODO: Estimate relative pose of two views. This can be subdivided into
//      - estimate the fundamental matrix F;
//      - compute the essential matrix E;
//      - recover rotation R and t.

// TODO: Reconstruct 3D points. The main task is
//      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

// TODO: Don't forget to
//          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
//          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
//            which can help you check if R and t are correct).
//       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
//       viewer will be notified to visualize the 3D points and update the view).
//       There are a few cases you should return 'false' instead, for example:
//          - function not implemented yet;
//          - input not valid (e.g., not enough points, point numbers don't match);
//          - encountered failure in any step.
