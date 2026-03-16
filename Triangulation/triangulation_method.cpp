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

    double u = norm_pts0[i].x();
    double v = norm_pts0[i].y();
    double up = norm_pts1[i].x();
    double vp = norm_pts1[i].y();


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

    Matrix F(3, 3, V.get_column(8).data());
    
    std::cout << "The F matrix: " << F << std::endl;

    Matrix UF(3, 3, 0.0);
    Matrix SF(3, 3, 0.0);
    Matrix VF(3, 3, 0.0);

    svd_decompose(F, UF, SF, VF);

    SF.set_row(2, {0.0, 0.0, 0.0    });

    Matrix F_constraint = UF * SF * VF.transpose();

    std::cout << "The F matrix after enforcing rank 2 constraint: " << F_constraint << std::endl;

    Matrix Final_F(3, 3, 0.0);
    Final_F = T1.transpose() * F_constraint * T0;

    std::cout << "The S matrix of F is: " << SF << std::endl;

    std::cout << "The final F matrix: " << Final_F << std::endl;

    Matrix P0Matrix(points_0.size(), 3, 1.0);
    Matrix P1Matrix(points_1.size(), 3, 1.0);

    for (size_t i = 0; i < points_0.size(); ++i) {
        P0Matrix.set_row(i, {points_0[i].x(), points_0[i].y(), 1.0});
        P1Matrix.set_row(i, {points_1[i].x(), points_1[i].y(), 1.0});
    }

    // std::cout << P1Matrix * Final_F * P0Matrix.transpose() << std::endl;

    Matrix33 F_unda = transpose(T1) * F_constraint * T0; //A 3×3 matrix encoding the epipolar geometry between the two views,
                                                    //It maps a point in image 0 to an epipolar line in image 1.
    Matrix33 K(fx, s,  cx,
                0, fy, cy,
                0,  0,  1);

    Matrix33 E = transpose(K) * F_unda* K;

    Matrix UE(3, 3, 0.0);
    Matrix SE(3, 3, 0.0);
    Matrix VE(3, 3, 0.0);

    svd_decompose(E, UE, SE, VE);

    Matrix33 W(0, -1.0, 0.0,
               1.0, 0.0, 0,
               0.0, 0, 1.0);

    Matrix R1 = UE * W * VE.transpose();
    Matrix R2 = UE * W.transpose() * VE.transpose();

    Vector t1 = UE.get_column(2);
    Vector t2 = -UE.get_column(2);
    

    std::cout << "The E matrix: " << E << std::endl;
    std::cout << "The first R matrix: " << R1 << std::endl;
    std::cout << "The alternative R matrix: " << R2 << std::endl;
    std::cout << "The positive t vector: " << t1 << std::endl;
    std::cout << "The negative t vector: " << t2 << std::endl;

    Matrix R_Used = R1;
    Vector t_Used = t2;

    Matrix M1 = K * Matrix34::identity();
    Matrix M2 = K * Matrix34(R_Used(0,0), R_Used(0,1), R_Used(0,2), t_Used[0],
                                R_Used(1,0), R_Used(1,1), R_Used(1,2), t_Used[1],
                                R_Used(2,0), R_Used(2,1), R_Used(2,2), t_Used[2]);
    std::cout << "The first M matrix: " << M1 << std::endl;
    std::cout << "The second M matrix: " << M2 << std::endl;


    Matrix A(4, 4, 0.0);
    
    for (int j = 0; j < points_0.size(); ++j) {
        Matrix A(4,4);
        A.set_row(0, points_0[j].x() * M1.get_row(2) - M1.get_row(0));
        A.set_row(1, points_0[j].y() * M1.get_row(2) - M1.get_row(1));
        A.set_row(2, points_1[j].x() * M2.get_row(2) - M2.get_row(0));
        A.set_row(3, points_1[j].y() * M2.get_row(2) - M2.get_row(1));
    
        Matrix U_A(4, 4, 0.0);
        Matrix S_A(4, 4, 0.0);
        Matrix V_A(4, 4, 0.0);
        svd_decompose(A, U_A, S_A, V_A);
        Vector4D X = V_A.get_column(3)/V_A.get_column(3)[3];
        points_3d.push_back(X.cartesian());
    }


    R = R_Used;
    t = t_Used;

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
