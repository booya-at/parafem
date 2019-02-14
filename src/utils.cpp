#include "utils.h"

namespace paraFEM {


inline bool is_finite(const Vector3 &x) {
    return ( (x - x).array() == (x - x).array()).all();
}

inline bool is_nan(const Vector3 &x){
    return !is_finite(x);
}


void check_nan(const Vector3 &vector, const char *message) {
    if (is_nan(vector)) {
        std::stringstream exception_message;
        exception_message << "Exception (" << vector.array() << message;
        std::string exception_string = exception_message.str();

        throw FemException(exception_string.c_str());
    }
}

}