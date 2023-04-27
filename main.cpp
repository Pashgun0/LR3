#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>

double rad(double ang)
{
    return ang * M_PI / 180;
}

struct Point
{
    double N;
    double E;

};

struct Dist
{
    std::string from;
    std::string to;
    std::string method;
    double value;
    int a;
    int b;
};

double getDist(Point& p1, Point& p2)
{
    return sqrt(pow(p1.N - p2.N, 2) + pow(p1.E - p2.E, 2));

}

double ctg(double angle)
{
    return 1 / tan(angle);
}

int main()
{
 
   
    Eigen::VectorXd dir(12);
    dir <<
    	0.00000000,		
    	23.173333333,	
    	50.108555556,	
    	0.00000000,		
    	48.397222222,	
    	132.946555556,	
    	0.00000000,		
    	48.520250000,	
    	72.399222222,	
    	0.00000000,		
    	45.340888889,	
    	104.543916667	
        ;
    Eigen::VectorXd dist(4);
    dist <<
        1305.530,
        778.820,
        701.850,
        1010.470;
    Eigen::MatrixXd S(8, 12);
    S <<
   -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0,
    0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    Eigen::VectorXd ang(8);
    ang = S * dir;

    Eigen::VectorXd b(8);
    b = ang * M_PI / 180;

    std::cout << ang << std::endl;
    Eigen::MatrixXd Qm = Eigen::MatrixXd::Zero(12, 12);
    for (int i = 0; i < 12; ++i)
    {
        Qm(i, i) = pow(1.4142, 2);

    }
    std::cout << "Qm =\n" << Qm << std::endl;

    Eigen::MatrixXd Qang = S * Qm * S.transpose();

    std::cout << "Qang =\n" << Qang << std::endl;

    Dist d0 = {"A", "C", "HD", 1305.530, 2, 3};
    Dist d1 = {"B", "D", "HD", 778.820, 2, 3};
    Dist d2 = {"B", "C", "HD", 701.850, 2, 3};
    Dist d3 = {"A", "D", "HD", 1010.470, 2, 3};
    std::vector<Dist> all_dist = {d0, d1, d2, d3};

    Eigen::MatrixXd Qs = Eigen::MatrixXd::Zero(dist.rows(), dist.rows());

    for (size_t i = 0; i < dist.rows(); i++)
    {
        double aa = all_dist.at(i).a;
        double bb = all_dist.at(i).b;
        double m = aa + bb * dist(i) / 1000.0;
        Qs(i,i) = pow(m /1000, 2); 
    }
    std::cout << "Qs = \n" << Qs <<std::endl;

    int N = 8 +4;
    int k = 4;
    int r = N - k;
    Eigen::VectorXd W(r);
    W(0) = ang(0) + ang(1) + ang(2) + ang(3) - 180;
    W(1) = ang(2) + ang(3) + ang(4) + ang(5) - 180;
    W(2) = ang(4) + ang(5) + ang(6) + ang(7) - 180;
    double num = sin(rad(ang(1))) *  sin(rad(ang(3))) * sin(rad(ang(5))) * sin(rad(ang(7)));
    double denum = sin(rad(ang(0))) *  sin(rad(ang(2))) * sin(rad(ang(4))) * sin(rad(ang(6)));
    W(3) = (1 - (num)/(denum)) * 180 / M_PI;

    Point pA{1077.000, 146.000};
    Point pB{1581.000, 663.000};
    double sAB = getDist(pA, pB);

    W(4) = 1 - sAB/dist(0) * sin(rad(ang(1)) + rad(ang(2)))/ sin(rad(ang(3)));
    W(5) = 1-  sAB/dist(1) * sin(rad(ang(0) + ang(7))) / sin(rad(ang(6)));
    W(6) = 1-  sAB/dist(2) * sin(rad(ang(0))) / sin(rad(ang(3)));
    W(7) = 1-  sAB/dist(3) * sin(rad(ang(1))) / sin(rad(ang(6)));

    std::cout << "ang=\n" << ang << std::endl;
    std::cout << "W =\n" << W << std::endl;

    for (size_t i = 0; i < 4; i++)
    {
        W(i) *= 3600;

    }


    for (size_t i = 4; i < W.rows(); i++)
    {
        W(i) *= 180 * 3600 / M_PI;
    }

    std::cout << "sAB =\n" << sAB << std::endl;
    std::cout << "W =\n" << W << std::endl;

    Eigen::MatrixXd B(8, 12);
    B <<
        1,                 1,                 1,                 1,             0,          0,              0,          0,                  0,                          0,                          0,                              0,                                 
        0,                 0,                 1,                 1,             1,          1,              0,          0,                  0,                          0,                          0,                              0,
        0,                 0,                 0,                 0,             1,          1,              1,          1,                  0,                          0,                          0,                              0,
        ctg(b(0)),         -1*ctg(b(1)),      ctg(b(2)),         -1*ctg(b(3)),  ctg(b(4)),  -1*ctg(b(5)),   ctg(b(6)),  -1*ctg(b(7)),       0,                          0,                          0,                              0,
        0,                 -1*ctg(b(1)+b(2)), -1*ctg(b(1)+b(2)), ctg(b(3)),     0,          0,              0,          0,                  1/dist(0)*180 *3600 / M_PI, 0,                          0,                              0,
        -1*ctg(b(0)+b(7)), 0,                 0,                 0,             0,          0,              ctg(b(6)),  -1*ctg(b(0)+b(7)),  0,                          1/dist(1)*180 *3600 / M_PI, 0,                              0,
        -1*ctg(b(0)),      0,                 0,                 ctg(b(3)),     0,          0,              0,          0,                  0,                          0,                           1/dist(2)*180 *3600 / M_PI,    0,
        0,                 -1*ctg(b(1)),      0,                 0,             0,          0,              ctg(b(6)),  0,                  0,                          0,                          0,                               1/dist(3)*180 *3600 / M_PI;
    
    std::cout << "B = \n" << B << std::endl;

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(12, 12);
    for (size_t i = 0; i < 8; i++)
    {
        for (size_t j = 0; j < 8; j++)
        {
            Q(i,j) = Qang(i,j);
        }
    }    

    for (size_t i = 8; i < 12; i++)
    {
        Q(i, i) = Qs(i-8,i-8);
    } 

    std::cout << "Q = \n" << Q << std::endl;
    Eigen::MatrixXd R = B * Q * B.transpose();
    std::cout << std::fixed << std::setprecision(9);
    std::cout << "R.det = \n" << R.determinant() << std::endl;
    Eigen::VectorXd K = -1 * R.inverse() * W;
    Eigen::VectorXd V = Q * B.transpose() * K;
    std::cout << "R = \n" << R << std::endl;
    std::cout << "K = \n" << K << std::endl;
    std::cout << "V = \n" << V << std::endl;
    std::cout << "B * V + W = \n" << B * V + W << std::endl; 

}