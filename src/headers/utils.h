#ifndef UTILS_H
#define UTILS_H

#include "Eigen/Core"

namespace parafem {

typedef Eigen::Vector3d Vector3;
typedef Eigen::Vector2d Vector2;

class FemException : public std::exception {
public:
    explicit FemException(const char * m) : message{m} {}
    virtual const char * what() const noexcept override {return message.c_str();}
private:
    std::string message = "";
};

bool is_finite(const Vector3 &x) ;

bool is_nan(const Vector3 &x);

void check_nan(const Vector3 &vector, const char *message);

}
#endif