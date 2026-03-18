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

double det3x3(const Matrix& M) {
    return 
        M(0,0) * (M(1,1)*M(2,2) - M(1,2)*M(2,1)) -
        M(0,1) * (M(1,0)*M(2,2) - M(1,2)*M(2,0)) +
        M(0,2) * (M(1,0)*M(2,1) - M(1,1)*M(2,0));
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

    Matrix U(norm_pts0.size(), norm_pts0.size(), 0.0);
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

    std::cout << "The S matrix of F is: " << SF << std::endl;


    Matrix P0Matrix(points_0.size(), 3, 1.0);
    Matrix P1Matrix(points_1.size(), 3, 1.0);

    for (size_t i = 0; i < points_0.size(); ++i) {
        P0Matrix.set_row(i, {points_0[i].x(), points_0[i].y(), 1.0});
        P1Matrix.set_row(i, {points_1[i].x(), points_1[i].y(), 1.0});
    }

    // std::cout << P1Matrix * Final_F * P0Matrix.transpose() << std::endl;

    Matrix33 F_unda = transpose(T1) * F_constraint * T0; 
    std::cout << "The final F matrix: " << F_unda << std::endl;
                                                    
    Matrix33 K(fx, s,  cx,
                0, fy, cy,
                0,  0,  1);

    Matrix33 E = transpose(K) * F_unda* K;

    Matrix UE(3, 3, 0.0);
    Matrix SE(3, 3, 0.0);
    Matrix VE(3, 3, 0.0);

    svd_decompose(E, UE, SE, VE);

    //  essential matrix constraint: It has rank 2 and is singular
    SE.set_row(0, {1, 0, 0});
    SE.set_row(1, {0, 1, 0});
    SE.set_row(2, {0, 0, 0});

    E = UE * SE * VE.transpose();

    Matrix33 W(0, -1.0, 0.0,
               1.0, 0.0, 0,
               0.0, 0, 1.0);

    Matrix R1 = UE * W * VE.transpose();
    Matrix R2 = UE * W.transpose() * VE.transpose();
    // Multiply by the determinant which is either +1 or -1
    R1 = det3x3(R1) * R1;
    R2 = det3x3(R2) * R2;

    Vector t1 = UE.get_column(2);
    Vector t2 = -UE.get_column(2);
    

    std::cout << "The E matrix: " << E << std::endl;
    std::cout << "The first R matrix: " << R1 << std::endl;
    std::cout << "The alternative R matrix: " << R2 << std::endl;
    std::cout << "The positive t vector: " << t1 << std::endl;
    std::cout << "The negative t vector: " << t2 << std::endl;

    std::vector<std::pair<Matrix, Vector3D>> candidates = {{R1,  t1}, {R1,  t2}, {R2,  t1}, {R2,  t2}
};
   
    int best_count = -1;  // even if points infront of cameras is 0, make it the best_count
    Matrix best_R;
    Vector3D best_t;
    std::vector<Vector3D> best_points;

    for (int i = 0; i < candidates.size(); i++) {

        Matrix R_test = candidates[i].first;
        Vector3D t_test = candidates[i].second;

        int count = 0;
        std::vector<Vector3D> temp_points;

        Matrix34 P1(1,0,0,0,
            0,1,0,0,
            0,0,1,0);

        Matrix34 P2(R_test(0,0), R_test(0,1), R_test(0,2), t_test[0],
            R_test(1,0), R_test(1,1), R_test(1,2), t_test[1],
            R_test(2,0), R_test(2,1), R_test(2,2), t_test[2]);

        Matrix M1 = K * P1;
        Matrix M2 = K * P2;
    
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
            Vector4D Xh = V_A.get_column(3)/V_A.get_column(3)[3];
       
            Vector3D P0 = Xh.cartesian();     // convert to 3D point
            temp_points.push_back(P0);       // appends 3D point to temp_points
            double z0 = P0.z();

            Vector3D P1 = R_test * P0 + t_test;  // transform point to camera 1's coordinate frame
            double z1 = P1.z();
            // check if point is infront of both cameras:
            if (z0 > 0 && z1 > 0) {
            count++;
            }
        }

    
    std::cout << "Combination " << i+1 << ": " << count << " points in front" << std::endl;
    
    // use the combination of R and t that gives the highest count of positive points
    if (count > best_count) {
        best_count = count;
        best_R = R_test;
        best_t = t_test;
        best_points = temp_points;
    }
}
    
    R = best_R;
    t = best_t;
    points_3d = best_points;

    std::cout << "The value for R is:" << R << "and the value for t is:" << t << std::endl;
    std::cout << "Best count: " << best_count << std::endl;

    return points_3d.size() > 0;
};

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
