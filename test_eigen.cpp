#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::VectorXd vec(3);
    vec << 1, 2, 3;
    std::cout << "Test vector: " << vec.transpose() << std::endl;
    return 0;
}