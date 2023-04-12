#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <cmath>

double rad(double ang)
{
    return ang * M_PI / 180;
}

struct Point
{
    double N;
    double E;

};

double getDist(Point& p1, Point& p2)
{
    return sqrt(pow(p1.N - p2.N, 2) + pow(p1.E - p2.E, 2));

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

    std::cout << ang << std::endl;
    Eigen::MatrixXd Qm = Eigen::MatrixXd::Zero(12, 12);
    for (int i = 0; i < 12; ++i)
    {
        Qm(i, i) = pow(1.4142, 2);

    }
    std::cout << "Qm =\n" << Qm << std::endl;

    Eigen::MatrixXd Qang = S * Qm * S.transpose();

    std::cout << "Qang =\n" << Qang << std::endl;

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


    std::cout << "sAB =\n" << sAB << std::endl;
    std::cout << "W =\n" << W * 3600 << std::endl;









}