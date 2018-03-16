#ifndef HLL_DEV_H__
#define HLL_DEV_H__

#include <random>
#include <algorithm>
#include <set>

namespace hll {

#ifdef ENABLE_HLL_DEVELOP
#pragma message("hll develop enabled (-DENABLE_HLL_DEVELOP)")
#else
namespace dev {
#endif

class hlldub_t: public hll_t {
    // hlldub_t inserts each value twice (forward and reverse)
    // and simply halves cardinality estimates.
public:
    template<typename... Args>
    hlldub_t(Args &&...args): hll_t(std::forward<Args>(args)...) {}
    INLINE void add(uint64_t hashval) {
        hll_t::add(hashval);
#ifndef NOT_THREADSAFE
        for(const uint32_t index(hashval & ((m()) - 1)), lzt(ffs(((hashval >> 1)|UINT64_C(0x8000000000000000)) >> (p() - 1)));
            core_[index] < lzt;
            __sync_bool_compare_and_swap(core_.data() + index, core_[index], lzt));
#else
        const uint32_t index(hashval & (m() - 1)), lzt(ffs(((hashval >> 1)|UINT64_C(0x8000000000000000)) >> (p() - 1)));
        core_[index] = std::min(core_[index], lzt);
#endif
    }
    double report() {
        sum();
        return this->creport();
    }
    double creport() {
        return hll_t::creport() * 0.5;
    }
    bool may_contain(uint64_t hashval) {
        return hll_t::may_contain(hashval) && core_[hashval & ((m()) - 1)] >= ffs(hashval >> p());
    }

    INLINE void addh(uint64_t element) {add(wang_hash(element));}
};

class dhll_t: public hll_t {
    // dhll_t is a bidirectional hll sketch which does not currently support set operations
    // It is based on the idea that the properties of a hll sketch work for both leading and trailing zeros and uses them as independent samples.
    std::vector<uint8_t, Allocator> dcore_;
public:
    template<typename... Args>
    dhll_t(Args &&...args): hll_t(std::forward<Args>(args)...),
                            dcore_(1ull << hll_t::p()) {
    }
    void sum() {
        uint64_t fcounts[64]{0};
        uint64_t rcounts[64]{0};
        const auto &core(hll_t::core());
        for(size_t i(0); i < core.size(); ++i) {
            // I don't this can be unrolled and LUT'd.
            ++fcounts[core[i]]; ++rcounts[dcore_[i]];
        }
        value_  = detail::calculate_estimate(fcounts, use_ertl_, m(), np_, alpha());
        value_ += detail::calculate_estimate(rcounts, use_ertl_, m(), np_, alpha());
        value_ *= 0.5;
        is_calculated_ = 1;
    }
    void add(uint64_t hashval) {
        hll_t::add(hashval);
#ifndef NOT_THREADSAFE
        for(const uint32_t index(hashval & ((m()) - 1)), lzt(ffs(((hashval >> 1)|UINT64_C(0x8000000000000000)) >> (p() - 1)));
            dcore_[index] < lzt;
            __sync_bool_compare_and_swap(dcore_.data() + index, dcore_[index], lzt));
#else
        const uint32_t index(hashval & (m() - 1)), lzt(ffs(((hashval >> 1)|UINT64_C(0x8000000000000000)) >> (p() - 1)));
        dcore_[index] = std::min(dcore_[index], lzt);
#endif
    }
    void addh(uint64_t element) {add(wang_hash(element));}
    bool may_contain(uint64_t hashval) {
        return hll_t::may_contain(hashval) && dcore_[hashval & ((m()) - 1)] >= ffs(((hashval >> 1)|UINT64_C(0x8000000000000000)) >> (p() - 1));
    }
};

static inline uint64_t finalize(uint64_t h) {
    // Murmurhash3 finalizer, for multiplying hash functions for seedhll_t and hllfilter_t.
    h ^= h >> 33;
    h *= 0xff51afd7ed558ccd;
    h ^= h >> 33;
    h *= 0xc4ceb9fe1a85ec53;
    h ^= h >> 33;
    return h;
}

class seedhll_t: public hll_t {
protected:
    uint64_t seed_; // 64-bit integers are xored with this value before passing it to a hash.
                          // This is almost free, in the content of 
public:
    template<typename... Args>
    seedhll_t(uint64_t seed, Args &&...args): hll_t(std::forward<Args>(args)...), seed_(seed) {
        if(seed_ == 0) LOG_WARNING("Note: seed is set to 0. No more than one of these at a time should have this value, and this is only for the purpose of multiplying hashes."
                                   " Also, if you are only using one of these at a time, don't use seedhll_t, just use hll_t and save yourself an xor per insertion"
                                   ", not to mention a 64-bit integer in space.");
    }
    seedhll_t(gzFile fp): hll_t() {
        this->read(fp);
    }
    void addh(uint64_t element) {
        element ^= seed_;
        add(wang_hash(element));
    }
    uint64_t seed() const {return seed_;}
    void write(const char *fn, bool write_gz) {
        if(write_gz) {
            gzFile fp = gzopen(fn, "wb");
            if(fp == nullptr) throw std::runtime_error("Could not open file.");
            this->write(fp);
            gzclose(fp);
        } else {
            std::FILE *fp = std::fopen(fn, "wb");
            if(fp == nullptr) throw std::runtime_error("Could not open file.");
            this->write(fileno(fp));
            std::fclose(fp);
        }
    }
    void write(gzFile fp) const {
        hll_t::write(fp);
        gzwrite(fp, &seed_, sizeof(seed_));
    }
    void read(gzFile fp) {
        hll_t::read(fp);
        gzread(fp, &seed_, sizeof(seed_));
    }
    void write(int fn) const {
        hll_t::write(fn);
        ::write(fn, &seed_, sizeof(seed_));
    }
    void read(int fn) {
        hll_t::read(fn);
        ::read(fn, &seed_, sizeof(seed_));
    }
    void read(const char *fn) {
        gzFile fp = gzopen(fn, "rb");
        this->read(fp);
        gzclose(fp);
    }

    template<typename T, typename Hasher=std::hash<T>>
    INLINE void adds(const T element, const Hasher &hasher) {
        static_assert(std::is_same_v<std::decay_t<decltype(hasher(element))>, uint64_t>, "Must return 64-bit hash");
        add(finalize(hasher(element) ^ seed_));
    }

#ifdef ENABLE_CLHASH
    template<typename Hasher=clhasher>
    INLINE void adds(const char *s, size_t len, const Hasher &hasher) {
        static_assert(std::is_same_v<std::decay_t<decltype(hasher(s, len))>, uint64_t>, "Must return 64-bit hash");
        add(finalize(hasher(s, len) ^ seed_));
    }
#endif
};

namespace detail {
inline std::set<uint64_t> seeds_from_seed(uint64_t seed, size_t size) {
    LOG_DEBUG("Initializing a vector of seeds of size %zu with a seed-seed of %" PRIu64 "\n", size, seed);
    std::mt19937_64 mt(seed);
    std::set<uint64_t> rset;
    while(rset.size() < size) rset.emplace(mt());
    return rset;
}
}

class hlf_t {
protected:
    // Consider templating this to extend to hlldub_ts as well.
    std::vector<seedhll_t> hlls_;
    mutable double value_;
    bool is_calculated_;
public:
    template<typename... Args>
    hlf_t(size_t size, uint64_t seedseed, Args &&... args): value_(0), is_calculated_(0) {
        auto sfs = detail::seeds_from_seed(seedseed, size);
        assert(sfs.size());
        hlls_.reserve(size);
        for(const auto seed: sfs)
            hlls_.emplace_back(seed, std::forward<Args>(args)...);
    }
    auto size() const {return hlls_.size();}
    auto m() const {return hlls_[0].size();}
    void write(const char *fn) const {
        gzFile fp = gzopen(fn, "wb");
        if(fp == nullptr) throw std::runtime_error("Could not open file.");
        this->write(fp);
        gzclose(fp);
    }
    void clear() {
        value_ = is_calculated_ = 0;
        for(auto &hll: hlls_) hll.clear();
    }
    void read(const char *fn) {
        gzFile fp = gzopen(fn, "rb");
        if(fp == nullptr) throw std::runtime_error("Could not open file.");
        this->read(fp);
        gzclose(fp);
    }
    void write(gzFile fp) const {
        uint64_t sz = hlls_.size();
        gzwrite(fp, &sz, sizeof(sz));
        for(const auto &hll: hlls_) {
            hll.write(fp);
        }
        gzclose(fp);
    }
    void read(gzFile fp) {
        uint64_t size;
        gzread(fp, &size, sizeof(size));
        hlls_.clear();
        while(hlls_.size() < size) hlls_.emplace_back(fp);
        gzclose(fp);
    }

    // This only works for hlls using 64-bit integers.
    // Looking ahead,consider templating so that the double version might be helpful.

    auto may_contain(uint64_t element) const {
#pragma message("Note: may_contain only works for the HyperLogFilter in the case of 64-bit integer insertions. One must hash a string to a 64-bit integer first in order to use it for this purpose.")
        return std::accumulate(hlls_.begin() + 1, hlls_.end(), hlls_.front().may_contain(wang_hash(element ^ hlls_.front().seed())),
                               [element](auto a, auto b) {
            return a && b.may_contain(wang_hash(element ^ b.seed()));
        });
    }
    void add(uint64_t val) {
        for(auto &hll: hlls_) hll.addh(val);
    }
    double creport() const {
        if(is_calculated_) return value_;
        double ret(hlls_[0].creport());
        for(size_t i(1); i < size(); ret += hlls_[i++].creport());
        ret /= static_cast<double>(size());
        value_ = ret;
        return value_;
    }
    double report() noexcept {
        if(is_calculated_) return value_;
        if(!hlls_[0].is_ready()) hlls_[0].sum();
        double ret(hlls_[0].report());
        for(size_t i(1); i < size(); ++i) {
            if(!hlls_[i].is_ready()) hlls_[i].sum();
            ret += hlls_[i].report();
        }
        ret /= static_cast<double>(size());
        return value_ = ret;
    }
    double med_report() noexcept {
        std::vector<double> values;
        values.reserve(size());
        for(auto hll: hlls_) values.emplace_back(hll.report());
        std::sort(values.begin(), values.end());
        if(size() & 1) return values[size() >> 1];
        return 0.5 * (values[size() >> 1] + values[(size() >> 1) - 1]);
    }
};

#if 0
    auto p = hll1.p();
    auto q = hll1.q();
    const double cAX = hll1.report();
    const double cBX = hll2.report();
    const double cABX = union_size(hll1, hll2);
    uint64_t AXBhalfsum = total_sum_of_stuff;
    uint64_t BXAhalfsum = total_sum_of_stuff;
    uint64_t cAXBhalf[64]{0}; // 
    uint64_t cBXAhalf[64]{0}; // 
    for(unsigned _q = 0; _q < q; ++_q)
        cAXBhalf[q] = num_larger1(q) + num_same(q) + num_larger2(q + 1);
        cAXBhalf[q] = num_larger2(q) + num_same(q) + num_larger1(q + 1);
        halfsums -= respective_halves[q];
    cAXBhalf[q] = AXBhalfsum;
    cBXAhalf[q] = BXAhalfsum;
    new_est_type = lambda p, q, data: estimate with p, but this time with q - 1
    cAXBhalf = new_est_type(countsAXBhalf, p, q) // And the omplment
    cA = cABX - cBX;
    cB = cABX - cAX;
    cX1 = (1.5 * cBX + 1.5*xAX - cBXAhalf - cAXBhalf);
    cX2 = 2.*(cBXAhalf + cAXBhalf) - 3.*cABX;
    return std::max(0, 0.5 * (cX1 + cX2));
/*
 *  Now go through and make 3 arrays:
 *
 */

    std::vector<int> countsAXBhalf(jointStatistic.getQ() + 1);
    std::vector<int> countsBXAhalf(jointStatistic.getQ() + 1);
    int sumAXBhalf= jointStatistic.getNumRegisters();
    int sumBXAhalf= jointStatistic.getNumRegisters();
    for (int q = 0; q < jointStatistic.getQ(); ++q) {
        countsAXBhalf[q] = jointStatistic.getLarger1Count(q) + jointStatistic.getEqualCount(q) + jointStatistic.getLarger2Count(q+1);
        countsBXAhalf[q] = jointStatistic.getLarger2Count(q) + jointStatistic.getEqualCount(q) + jointStatistic.getLarger1Count(q+1);
        sumAXBhalf -= countsAXBhalf[q];
        sumBXAhalf -= countsBXAhalf[q];
    }
    countsAXBhalf[jointStatistic.getQ()] = sumAXBhalf;
    countsBXAhalf[jointStatistic.getQ()] = sumBXAhalf;
    const MaxLikelihoodEstimator estimator2(jointStatistic.getP(), jointStatistic.getQ()-1);

    const double cardinalityAXBhalf = estimator2(countsAXBhalf);
    const double cardinalityBXAhalf = estimator2(countsBXAhalf);

    cardinalityA = cardinalityABX - cardinalityBX;
    cardinalityB = cardinalityABX - cardinalityAX;
    double cardinalityX1 = 1.5*cardinalityBX + 1.5*cardinalityAX - cardinalityBXAhalf - cardinalityAXBhalf;
    double cardinalityX2 = 2.*(cardinalityBXAhalf + cardinalityAXBhalf) - 3*cardinalityABX;

cardinalityX = std::max(0., 0.5*(cardinalityX1 + cardinalityX2));
#endif

#ifdef ENABLE_HLL_DEVELOP
#pragma message("hll develop enabled (-DENABLE_HLL_DEVELOP)")
#else
} // namespace dev
#endif
}

#endif // #ifndef HLL_DEV_H__
