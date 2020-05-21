#include "BigInt.h"
#include "math.h"
#include <iostream>

void test_pow() {
    BigInt x = 10;
    BigInt y = 20;
    BigInt z = power(x, y);

    std::cout << x.ToString() << " ** " << y.ToString() << " = " << z.ToString() << std::endl;
    //fprintf(stderr, z.ToString());
    fprintf(stderr, "test result = %s, and %d\n", z.ToString(), pow(10, 20));
    BigInt a = 2;
    BigInt b = 5;
    BigInt c = power(a, b);
    std::cout << a.ToString() << " ** " << b.ToString() << " = " << c.ToString() << std::endl;
    std::cout << "5**5 = " << power(5, 5) << std::endl;
}

int main() {
	test_pow();
}
