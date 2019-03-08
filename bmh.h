#ifndef SKETCH_BAG_MINHASH__
#define SKETCH_BAG_MINHASH__
#include "common.h"
#include "bagminhash/c++/weighted_minwise_hashing.hpp"
#include "bbmh.h"

namespace sketch {

namespace minhash {
// Include here to wrap in a namespace, since there's no namespacing in bagminhash
using FloatBagMinHasher=BagMinHasher2<FloatWeightDiscretization, XXHash64, common::Allocator>;
using BinaryBagMinHasher=BagMinHasher2<BinaryWeightDiscretization, XXHash64, common::Allocator>;

}

}

#endif /* SKETCH_BAG_MINHASH__ */
