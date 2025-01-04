#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <chrono>
#include <random>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <complex>
#include <cstdlib>
#include <random>
#include <numeric>

const int degree = 8;
typedef std::complex<double> dcomp;
typedef std::complex<double> Complex;
typedef std::vector<Complex> Poly;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(-10.0, 10.0);

std::pair<double, double> computeVariance(const std::vector<double> &data)
{
    std::pair<bool, int> result;
    int n = data.size();
    // if (n <= 1) return 0.0; // 方差至少需要两个数据点

    // 计算平均值
    double sum = std::accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / n;

    // 计算方差
    double variance = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double diff = data[i] - mean;
        variance += diff * diff;
    }
    variance /= n;
    result.first = mean;
    result.second = variance;
    std::cout << "Mean: " << mean << ", Variance: " << variance << std::endl;
    return result;
}

Complex evaluate(const Poly &poly, const Complex &x)
{
    Complex result = 0;
    for (auto it = poly.rbegin(); it != poly.rend(); ++it)
    {
        result = result * x + *it;
    }
    return result;
}

// 计算多项式的导数
Poly derivative(const Poly &poly)
{
    Poly deriv(poly.size() - 1);
    for (size_t i = 1; i < poly.size(); ++i)
    {
        deriv[i - 1] = poly[i] * Complex(i, 0);
    }
    return deriv;
}

// Aberth方法求解多项式根
std::vector<Complex> aberth(const Poly &poly)
{
    size_t n = poly.size() - 1;
    Poly deriv = derivative(poly);
    std::vector<Complex> roots(n, 0);

    // 初始化根的初始猜测值s
    const double PI = 3.14159265358979323846;
    for (size_t i = 0; i < n; ++i)
    {
        roots[i] = std::polar(1.0, 2 * PI * i / n);
    }

    std::vector<Complex> delta(n);
    const double tolerance = 1e-5;
    bool converged = false;

    while (!converged)
    {
        converged = true;
        for (size_t i = 0; i < n; ++i)
        {
            Complex f_val = evaluate(poly, roots[i]);
            Complex f_prime_val = evaluate(deriv, roots[i]);
            Complex sum = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    sum += 1.0 / (roots[i] - roots[j]);
                }
            }

            delta[i] = f_val / (f_prime_val - f_val * sum);

            if (std::abs(delta[i]) > tolerance)
            {
                converged = false;
            }
        }

        for (size_t i = 0; i < n; ++i)
        {
            roots[i] -= delta[i];
        }
    }

    return roots;
}

dcomp poly(double A[], int n, dcomp x)
{
    dcomp y = pow(x, n);
    for (int i = 0; i < n; i++)
        y += A[i] * pow(x, (n - i - 1));
    return y;
}

dcomp *polyroot(double A[], int n, double tolerance)
{
    int iterations = 100;
    dcomp z = dcomp(0.4, 0.9);
    int size = sizeof(z);
    dcomp *R;
    dcomp *R_pre;
    R = (dcomp *)malloc(size * n);
    for (int i = 0; i < n; i++)
        R[i] = pow(z, i);
    for (int i = 0; i < iterations; i++)
    {
        R_pre = R;
        for (int j = 0; j < n; j++)
        {
            dcomp B = poly(A, n, R[j]);
            for (int k = 0; k < n; k++)
            {
                if (k != j)
                    B /= R[j] - R[k];
            }
            R[j] -= B;
        }
        for (int j = 0; j < n; j++)
        {
            if (abs((R[j] - R_pre[j])) < tolerance)
                break;
        }
    }
    return R;
}

int Fac(int n)
{
    if (n == 0 || n == 1)
    {
        return 1;
    }
    else
    {
        return Fac(n - 1) * n;
    }
}

// 定义多项式类
class Polynomial
{
public:
    Polynomial(int degree, std::vector<double> coefficients) : degree_(degree), coefficients_(coefficients) {}

    double operator()(double x) const
    {
        double result = 0;
        for (int i = 0; i <= degree_; ++i)
        {
            result += coefficients_[i] * std::pow(x, i);
        }
        return result;
    }

    void print() const
    {
        std::stringstream ss;
        ss << coefficients_[0];
        for (int i = 1; i <= degree_; ++i)
        {
            if (coefficients_[i] >= 0)
            {
                ss << " + " << coefficients_[i] << "x^" << i;
            }
            else
            {
                ss << " - " << -coefficients_[i] << "x^" << i;
            }
        }
        std::cout << "Polynomial : " << ss.str() << std::endl;
    }

    double derivative(double x) const
    {
        double result = 0;
        for (int i = 1; i <= degree_; ++i)
        {
            result += i * coefficients_[i] * std::pow(x, i - 1);
        }
        return result;
    }

    double derivative(double x, int n) const
    {
        double result = 0;
        for (int i = n; i <= degree_; ++i)
        {
            result += Fac(i) * coefficients_[i] / Fac(i - n) * std::pow(x, i - n);
        }
        return result;
    }

    double con_val(double x) const
    {
        int n = 1;
        double result0 = 0;
        for (int i = n; i <= degree_; ++i)
        {
            result0 += fabs(Fac(i) * coefficients_[i] / Fac(i - n) * std::pow(x, i - n));
        }
        double result1 = 0;
        for (int i = 1; i <= degree_; ++i)
        {
            result1 += fabs(i * coefficients_[i] * std::pow(x, i - 1));
        }
        double result = result1 / result0;
        return result;
    }

    double fd(double x) const
    {
        return derivative(x, 0) / derivative(x, 1);
    }

    float muller(float a, float b, float c)
    {
        const int MAX = 100;
        int i;
        float result;
        for (i = 0;; ++i)
        {
            // Calculating various constants required
            // to calculate x3
            float f1 = derivative(a, 0);
            float f2 = derivative(a, 0);
            float f3 = derivative(a, 0);
            float d1 = f1 - f3;
            float d2 = f2 - f3;
            float h1 = a - c;
            float h2 = b - c;
            float a0 = f3;
            float a1 = (((d2 * pow(h1, 2)) - (d1 * pow(h2, 2))) / ((h1 * h2) * (h1 - h2)));
            float a2 = (((d1 * h2) - (d2 * h1)) / ((h1 * h2) * (h1 - h2)));
            float x = ((-2 * a0) / (a1 + fabs(sqrt(a1 * a1 - 4 * a0 * a2))));
            float y = ((-2 * a0) / (a1 - fabs(sqrt(a1 * a1 - 4 * a0 * a2))));
            // Taking the root which is closer to x2
            if (x >= y)
                result = x + c;
            else
                result = y + c;
            // checking for resemblance of x3 with x2 till
            // two decimal places
            float m = result * 100;
            float n = c * 100;
            m = floor(m);
            n = floor(n);
            if (m == n)
                break;
            a = b;
            b = c;
            c = result;
            if (i > MAX)
            {
                printf("Root can't be found using Muller method");
                break;
            }
        }
        if (i <= MAX)
            return result;
    }

    double po(double a, double b)
    {
        double eps = 0.00001;
        double x;
        int i = 0;
        do
        {
            x = (a * derivative(b, 0) - b * derivative(a, 0)) / (derivative(b, 0) - derivative(a, 0));
            if (derivative(x, 0) * derivative(a, 0) > 0)
                a = x;
            else
                b = x;
            // std::cout << x << ',' << derivative(x, 0) << std::endl;
            // std::cout << a << ',' << b << std::endl;
            i++;
            if (i == 100)
                break;
        } while (fabs(derivative(x, 0)) > eps);
        std::cout << "迭代次数为：" << i << std::endl;
        std::cout << "x : " << x << std::endl;
        std::cout << "fx : " << derivative(x, 0) << std::endl;
        return x;
    }

    float dichotomy(float a, float b, float e)
    {
        if (derivative(a, 0) * derivative(b, 0) > 0)
        {
            // 查找算法失效
            std::cout << "dichotomy失败" << std::endl;
            return 0.0f;
        }
        int i = 0;
        while (true)
        {
            // 计算中间值
            float x = 0.5f * (a + b);
            // 获取函数值的绝对值
            float f = fabs(derivative(x, 0));
            float tol = fabs(b - a);
            // 判断是否跳出循环
            if (f < e || tol < e)
            {
                std::cout << "迭代次数为：" << i << std::endl;
                std::cout << "x : " << x << std::endl;
                std::cout << "fx : " << derivative(x, 0) << std::endl;
                return x;
            }
            // 开始折半搜索查找
            if (derivative(x, 0) * derivative(a, 0) < 0)
            {
                b = x;
            }
            else
            {
                a = x;
            }
            i++;
        }
    }

    double newton(double x0, double e)
    {
        double a = x0;
        double x = x0 - fd(x0);
        int i = 0;
        // double judge = double(a) - double(x);
        // std::cout << judge << std::endl;
        // std::cout << fabs(judge) << std::endl;
        while (fabs(double(x) - double(a)) > e)
        {
            // std::cout << a << std::endl;
            a = x;
            i++;
            x = x - fd(x);
            if (i > 1000)
            {
                std::cout << "迭代失败！" << std::endl;
                return x;
            }
        }
        std::cout << "迭代次数为：" << i << std::endl;
        std::cout << "x : " << x << std::endl;
        std::cout << "fx : " << derivative(x, 0) << std::endl;
        return x;
    }

    int adj()
    {
        Eigen::Matrix<double, degree, degree> matrix_55;
        // 复数动态矩阵
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> matrix_eigenvalues;
        // 同样测试 12345
        for (int i = 0; i < degree; i++)
        {
            for (int j = 0; j < degree; j++)
            {
                if (j == degree - 1)
                {
                    matrix_55(i, j) = -coefficients_[i] / coefficients_[degree];
                }
                else if (i - j == 1)
                {
                    matrix_55(i, j) = 1;
                }
                else
                {
                    matrix_55(i, j) = 0;
                }
            }
        }
        // matrix_55 << 0, 0, 0, 0, -5,
        //     1, 0, 0, 0, -4,
        //     0, 1, 0, 0, -3,
        //     0, 0, 1, 0, -2,
        //     0, 0, 0, 1, -1;
        double condition_number = matrix_55.norm() * matrix_55.inverse().norm();

        std::cout << "矩阵的条件数为: " << condition_number << std::endl;
        std::cout << "matrix_55: " << std::endl
                  << matrix_55 << std::endl
                  << std::endl;

        matrix_eigenvalues = matrix_55.eigenvalues();
        // std::cout<<"matrix_55.eigenvalues: "<<std::endl<<matrix_55.eigenvalues()<<std::endl;
        std::cout << "matrix_eigenvalues: " << std::endl;
        std::cout << matrix_eigenvalues << std::endl;
        return 0;
    }

    int Durand_Kerner()
    {
        const double tolerance = 1e-10;
        double A[degree - 1];
        for (int a = 0; a < degree; a = a + 1)
        {
            A[degree - 1 - a] = coefficients_[a] / coefficients_[degree];
        }
        std::cout << "The coefficients of the polynomial are: ";
        for (int i = 0; i < degree; i = i + 1)
        {
            std::cout << A[degree - 1 - i] << " ";
        }
        std::cout << std::endl;
        dcomp *R = polyroot(A, degree, tolerance);
        for (int i = 0; i < degree; i++)
            std::cout << "R[" << i << "] = " << R[i] << std::endl;
        return 0;
    }

    int abe()
    {
        Poly poly;
        // for (int i = 0; i < degree + 1; i = i + 1)
        // {
        //     std::cout << coefficients_[i] << " ";
        // }
        for (int i = 0; i < degree + 1; ++i)
        {
            poly.push_back(coefficients_[i]);
        }
        // for (int i = 0; i < degree + 1; i = i + 1)
        // {
        //     std::cout << poly[i] << " ";
        // }

        std::vector<Complex> roots = aberth(poly);

        // 输出结果
        std::cout << "Roots: \n";
        for (const auto &root : roots)
        {
            std::cout << root << "\n";
        }

        return 0;
    }

private:
    int degree_;
    std::vector<double> coefficients_;
};

// 生成随机多项式
Polynomial generateRandomPolynomial(int degree)
{
    std::vector<double> coefficients(degree + 1);
    for (int i = 0; i < degree + 1; ++i)
    {
        // double r1 = circle1.Radius();
        // double r2 = circle2.Radius();
        // double xX = dirX.X();
        // double yX = dirX.Y();
        // double zX = dirX.Z();
        // double xY = dirY.X();
        // double yY = dirY.Y();
        // double zY = dirY.Z();
        // double x0 = center1.X();
        // double y0 = center1.Y();
        // double z0 = center1.Z();
        // double m = center1.Dot(circle1.Axes().DirectionX());
        // double n = center1.Dot(circle1.Axes().DirectionY());

        // double r1 = 2;
        // double r2 = 1;
        // double xX = 1, yX = 0, zX = 0;
        // double xY = 0, yY = std::pow(2, 1/2), zY = std::pow(2, 1/2);
        // double x0 = 245, y0 = 125, z0 = 786;
        // double r1 = 2;
        // double r2 = 1;
        // double xX = 1, yX = 0, zX = 0;
        // double xY = 0, yY = 0, zY = 1;
        // double x0 = 1, y0 = -1, z0 = 0;
        double r1 = 2;
        double r2 = 1;
        double xX = 1, yX = 0, zX = 0;
        double xY = 0, yY = 0, zY = 1;
        double x0 = 0, y0 = 1, z0 = 0;
        double m = x0 * xX + y0 * yX + z0 * zX;
        double n = x0 * xY + y0 * yY + z0 * zY;
        // matrix_eigenvalues:
        // (6.15854e-10,0.000135708)
        // (6.15854e-10,-0.000135708)
        //         (0.000135707,0)
        //         (-0.000135708,0)
        //                     (-1,0)
        //                     (1,0)
        //                     (-1,0)
        //                     (1,0)
        // adj38247.4微秒
        // The coefficients of the polynomial are: 0 -2 0 1 0 0 0 0
        // R[0] = (1,0)
        // R[1] = (6.85585e-18,-4.01048e-17)
        // R[2] = (-3.84678e-17,-6.45443e-18)
        // R[3] = (-1,-2.17108e-91)
        // R[4] = (-1,2.05166e-91)
        // R[5] = (1,-6.2168e-17)
        // R[6] = (3.32673e-17,2.56609e-18)
        // R[7] = (-2.55806e-18,3.18851e-17)
        // Durand_Kerner26043.9微秒

        // m^2*r1^2*(r1^2*xY^2 + r1^2*yY^2) - r2^2*r1^4*zX^2*zY^2
        double a0 = m * m * r1 * r1 * (r1 * r1 * xY * xY + r1 * r1 * yY * yY) - r2 * r2 * std::pow(r1, 4) * zX * zX * zY * zY;

        // m^2*r1^2*(2*r1^2*xX*xY + 2*r1^2*yX*yY) + 2*r2^2*r1^4*zX*zY^3 - 2*r2^2*r1^4*zX^3*zY - 2*m*n*r1^2*(r1^2*xY^2 + r1^2*yY^2)
        double b1 = m * m * r1 * r1 * (2 * r1 * r1 * xX * xY + 2 * r1 * r1 * yX * yY) + 2 * r2 * r2 * std::pow(r1, 4) * zX * std::pow(zY, 3) - 2 * r2 * r2 * std::pow(r1, 4) * std::pow(zX, 3) * zY - 2 * m * n * r1 * r1 * (r1 * r1 * xY * xY + r1 * r1 * yY * yY);

        // m^2*r1^2*(2*r1*x0*xY + 2*r1*y0*yY) - 2*r2^2*r1^3*z0*zX^2*zY
        double b0 = m * m * r1 * r1 * (2 * r1 * x0 * xY + 2 * r1 * y0 * yY) - 2 * r2 * r2 * std::pow(r1, 3) * z0 * zX * zX * zY;

        // m^2*r1^2*(r1^2*xX^2 + r1^2*yX^2) + n^2*r1^2*(r1^2*xY^2 + r1^2*yY^2) - r2^2*r1^4*zX^4 - r2^2*r1^4*zY^4
        //- 2*m*n*r1^2*(2*r1^2*xX*xY + 2*r1^2*yX*yY) + 4*r2^2*r1^4*zX^2*zY^2
        double c2 = m * m * r1 * r1 * (r1 * r1 * xX * xX + r1 * r1 * yX * yX) + n * n * r1 * r1 * (r1 * r1 * xY * xY + r1 * r1 * yY * yY) - r2 * r2 * std::pow(r1, 4) * std::pow(zX, 4) - r2 * r2 * std::pow(r1, 4) * std::pow(zY, 4) - 2 * m * n * r1 * r1 * (2 * r1 * r1 * xX * xY + 2 * r1 * r1 * yX * yY) + 4 * r2 * r2 * std::pow(r1, 4) * zX * zX * zY * zY;

        // m^2*r1^2*(2*r1*x0*xX + 2*r1*y0*yX) - 2*r2^2*r1^3*z0*zX^3 - 2*m*n*r1^2*(2*r1*x0*xY + 2*r1*y0*yY) + 4*r2^2*r1^3*z0*zX*zY^2
        double c1 = m * m * r1 * r1 * (2 * r1 * x0 * xX + 2 * r1 * y0 * yX) - 2 * r2 * r2 * std::pow(r1, 3) * z0 * std::pow(zX, 3) - 2 * m * n * r1 * r1 * (2 * r1 * x0 * xY + 2 * r1 * y0 * yY) + 4 * r2 * r2 * std::pow(r1, 3) * z0 * zX * zY * zY;

        // m^2*r1^2*(x0^2 + y0^2) - r2^2*r1^2*z0^2*zX^2
        double c0 = m * m * r1 * r1 * (x0 * x0 + y0 * y0) - r2 * r2 * r1 * r1 * z0 * z0 * zX * zX;

        // n^2*r1^2*(2*r1^2*xX*xY + 2*r1^2*yX*yY) - 2*r2^2*r1^4*zX*zY^3 + 2*r2^2*r1^4*zX^3*zY - 2*m*n*r1^2*(r1^2*xX^2 + r1^2*yX^2)
        double d3 = n * n * r1 * r1 * (2 * r1 * r1 * xX * xY + 2 * r1 * r1 * yX * yY) - 2 * r2 * r2 * std::pow(r1, 4) * zX * std::pow(zY, 3) + 2 * r2 * r2 * std::pow(r1, 4) * std::pow(zX, 3) * zY - 2 * m * n * r1 * r1 * (r1 * r1 * xX * xX + r1 * r1 * yX * yX);

        // n^2*r1^2*(2*r1*x0*xY + 2*r1*y0*yY) - 2*r2^2*r1^3*z0*zY^3 - 2*m*n*r1^2*(2*r1*x0*xX + 2*r1*y0*yX) + 4*r2^2*r1^3*z0*zX^2*zY
        double d2 = n * n * r1 * r1 * (2 * r1 * x0 * xY + 2 * r1 * y0 * yY) - 2 * r2 * r2 * std::pow(r1, 3) * z0 * std::pow(zY, 3) - 2 * m * n * r1 * r1 * (2 * r1 * x0 * xX + 2 * r1 * y0 * yX) + 4 * r2 * r2 * std::pow(r1, 3) * z0 * zX * zX * zY;

        // 2*r2^2*r1^2*z0^2*zX*zY - 2*m*n*r1^2*(x0^2 + y0^2)
        double d1 = 2 * r2 * r2 * r1 * r1 * z0 * z0 * zX * zY - 2 * m * n * r1 * r1 * (x0 * x0 + y0 * y0);

        // n^2*r1^2*(r1^2*xX^2 + r1^2*yX^2) - r2^2*r1^4*zX^2*zY^2
        double e4 = n * n * r1 * r1 * (r1 * r1 * xX * xX + r1 * r1 * yX * yX) - r2 * r2 * std::pow(r1, 4) * zX * zX * zY * zY;

        // n^2*r1^2*(2*r1*x0*xX + 2*r1*y0*yX) - 2*r2^2*r1^3*z0*zX*zY^2
        double e3 = n * n * r1 * r1 * (2 * r1 * x0 * xX + 2 * r1 * y0 * yX) - 2 * r2 * r2 * std::pow(r1, 3) * z0 * zX * zY * zY;

        // n^2*r1^2*(x0^2 + y0^2) - r2^2*r1^2*z0^2*zY^2
        double e2 = n * n * r1 * r1 * (x0 * x0 + y0 * y0) - r2 * r2 * r1 * r1 * z0 * z0 * zY * zY;

        coefficients[0] = (a0 + c0) * (a0 + c0) - b0 * b0;
        coefficients[1] = 2 * c1 * (a0 + c0) - 2 * b0 * (b1 + d1);
        coefficients[2] = 2 * b0 * (b0 - d2) - (b1 + d1) * (b1 + d1) - 2 * (a0 + c0) * (2 * a0 + c0 - c2 - e2) + b0 * b0 + c1 * c1;
        coefficients[3] = 2 * (b1 + d1) * (b0 - d2) - 2 * (a0 + c0) * (c1 - e3) + 2 * b0 * (b1 - d3) - 2 * c1 * (2 * a0 + c0 - c2 - e2) + 2 * b0 * (b1 + d1);
        coefficients[4] = 2 * (b1 + d1) * (b1 - d3) - 2 * b0 * (b0 - d2) - 2 * c1 * (c1 - e3) + (2 * a0 + c0 - c2 - e2) * (2 * a0 + c0 - c2 - e2) + 2 * (a0 + c0) * (a0 - c2 + e4) + (b1 + d1) * (b1 + d1) - (b0 - d2) * (b0 - d2);
        coefficients[5] = 2 * c1 * (a0 - c2 + e4) - 2 * b0 * (b1 - d3) - 2 * (b0 - d2) * (b1 - d3) - 2 * (b1 + d1) * (b0 - d2) + 2 * (c1 - e3) * (2 * a0 + c0 - c2 - e2);
        coefficients[6] = (b0 - d2) * (b0 - d2) - 2 * (a0 - c2 + e4) * (2 * a0 + c0 - c2 - e2) - 2 * b1 * (b1 - d3) - (b1 - d3) * (b1 - d3) + (c1 - e3) * (c1 - e3);
        coefficients[7] = 2 * (b0 - d2) * (b1 - d3) - 2 * (c1 - e3) * (a0 - c2 + e4);
        coefficients[8] = (a0 - c2 + e4) * (a0 - c2 + e4) + (b1 - d3) * (b1 - d3);

        // if (i == 0)
        // {
        //     coefficients[i] = 0.1 * (rand() % 20 - 10); // 生成-10到10之间的随机数
        // }
        // else
        // {
        //     coefficients[i] = i * (rand() % 20 - 10); // 生成-10到10之间的随机数
        // }
        // if (i == degree)
        // {
        //     coefficients[i] = (1 + i) * (1 + i) * (1 + i) * dis(gen); // 生成-10到10之间的随机数
        // }
        // else
        // {
        //     coefficients[i] = dis(gen); // 生成-10到10之间的随机数
        // }
        // coefficients[i] = rand() % 20 - 10; // 生成-10到10之间的随机数
    }
    return Polynomial(degree, coefficients);
}

int main()
{
    double x = 2.00001;
    std::vector<double> duration0;
    std::vector<double> duration1;
    std::vector<double> duration2;
    std::vector<double> duration3;
    std::vector<double> duration4;
    std::vector<double> duration5;
    std::vector<double> duration6;
    const int epoch = 100;

    for (int i = 0; i <= epoch; ++i)
    {
        Polynomial polynomial = generateRandomPolynomial(degree);
        polynomial.print(); // 打印多项式
        std::cout << "多项式的值为: " << polynomial(x) << std::endl;
        std::cout << "导的值为: " << polynomial.derivative(x, 0) << std::endl;

        auto beforeTime0 = std::chrono::steady_clock::now();
        double sol;
        sol = polynomial.newton(2.5, 0.00001);
        auto afterTime0 = std::chrono::steady_clock::now();
        double duration_microsecond = std::chrono::duration<double, std::micro>(afterTime0 - beforeTime0).count();
        std::cout << "newton" << duration_microsecond << "微秒" << std::endl;
        duration0.push_back(duration_microsecond);
        std::cout << "cn: " << polynomial.con_val(sol) << std::endl;

        auto beforeTime_1 = std::chrono::steady_clock::now();
        polynomial.muller(-10.0, 5.0, 10.0);
        auto afterTime_1 = std::chrono::steady_clock::now();
        double duration_microsecond1 = std::chrono::duration<double, std::micro>(afterTime_1 - beforeTime_1).count();
        std::cout << "muller" << duration_microsecond1 << "微秒" << std::endl;
        duration1.push_back(duration_microsecond1);

        auto beforeTime_2 = std::chrono::steady_clock::now();
        polynomial.po(-20.0, 20.0);
        auto afterTime_2 = std::chrono::steady_clock::now();
        double duration_microsecond2 = std::chrono::duration<double, std::micro>(afterTime_2 - beforeTime_2).count();
        std::cout << "po" << duration_microsecond2 << "微秒" << std::endl;
        duration2.push_back(duration_microsecond2);

        auto beforeTime_3 = std::chrono::steady_clock::now();
        polynomial.dichotomy(-20.0, 20.0, 0.00001);
        auto afterTime_3 = std::chrono::steady_clock::now();
        double duration_microsecond3 = std::chrono::duration<double, std::micro>(afterTime_3 - beforeTime_3).count();
        std::cout << "dichotomy" << duration_microsecond3 << "微秒" << std::endl;
        duration3.push_back(duration_microsecond3);

        auto beforeTime_4 = std::chrono::steady_clock::now();
        polynomial.adj();
        auto afterTime_4 = std::chrono::steady_clock::now();
        double duration_microsecond4 = std::chrono::duration<double, std::micro>(afterTime_4 - beforeTime_4).count();
        std::cout << "adj" << duration_microsecond4 << "微秒" << std::endl;
        duration4.push_back(duration_microsecond4);

        auto beforeTime_5 = std::chrono::steady_clock::now();
        polynomial.Durand_Kerner();
        auto afterTime_5 = std::chrono::steady_clock::now();
        double duration_microsecond5 = std::chrono::duration<double, std::micro>(afterTime_5 - beforeTime_5).count();
        std::cout << "Durand_Kerner" << duration_microsecond5 << "微秒" << std::endl;
        duration5.push_back(duration_microsecond5);

        auto beforeTime_6 = std::chrono::steady_clock::now();
        polynomial.abe();
        auto afterTime_6 = std::chrono::steady_clock::now();
        double duration_microsecond6 = std::chrono::duration<double, std::micro>(afterTime_6 - beforeTime_6).count();
        std::cout << "aberth" << duration_microsecond6 << "微秒" << std::endl;
        duration6.push_back(duration_microsecond6);
    }

    computeVariance(duration0);
    computeVariance(duration1);
    computeVariance(duration2);
    computeVariance(duration3);
    computeVariance(duration4);
    computeVariance(duration5);
    computeVariance(duration6);

    return 0;
}
