#include "bmh.h"
using namespace sketch;
using namespace mh;

int main() {
    FloatBagMinHasher fbhm(2000);
    for(size_t i = 100000; --i;fbhm.addh((uint64_t(rand()) << 32) | rand()));
    auto f = FinalDivBBitMinHash(fbhm.finalize());
}
