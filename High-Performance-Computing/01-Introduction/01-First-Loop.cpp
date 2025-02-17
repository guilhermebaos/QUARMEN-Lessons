#include <iostream>
#include <math.h>

int main () {

    for (unsigned i = 0; i < 5; i++) {
        double x = M_PI * i;
        std::cout << i << " -> " << x << std::endl;
    }

    std::cout << std::endl;

    unsigned j = 0;
    unsigned k = 0;
    while (j < 4) {
        std::cout << "j = " << ++j << "   " << "k = " << k++ << std::endl;

        if (k < 2) {
            std::cout << "We have k < 2" << std::endl;
        } else {
            std::cout << "We have k >= 2" << std::endl;
        }
    }
}