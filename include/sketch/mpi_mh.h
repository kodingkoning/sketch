#pragma once
#include <mutex>
//#include <queue>
#include "hll.h" // For common.h and clz functions
#include "fixed_vector.h"
#include <unordered_map>
#include "isz.h"


/*
 * TODO: support minhash using sketch size and a variable number of hashes.
 * Implementations: KMinHash (multiple hash functions)
 *                  RangeMinHash (the lowest range in a set, but only uses one hash)
 *                  HyperMinHash
 */

namespace sketch {
inline namespace mpi_minhash {

#define SET_SKETCH(x) const auto &sketch() const {return x;} auto &sketch() {return x;}

namespace detail {



static constexpr double HMH_C = 0.169919487159739093975315012348;

template<typename FType, typename=typename std::enable_if<std::is_floating_point<FType>::value>::type>
inline FType beta(FType v) {
    FType ret = -0.370393911 * v;
    v = std::log1p(v);
    FType poly = v * v;
    ret += 0.070471823 * v;
    ret += 0.17393686 * poly;
    poly *= v;
    ret += 0.16339839 * poly;
    poly *= v;
    ret -= 0.09237745 * poly;
    poly *= v;
    ret += 0.03738027 * poly;
    poly *= v;
    ret += -0.005384159 * poly;
    poly *= v;
    ret += 0.00042419 * poly;
    return ret;
}

} // namespace detail

template<typename T, typename Cmp=std::greater<T>>
class AbstractMinHash {
protected:
    uint64_t ss_;
    AbstractMinHash(uint64_t sketch_size): ss_(sketch_size) {}
    auto &sketch();
    const auto &sketch() const;
    auto nvals() const {return ss_;}
    std::vector<T, common::Allocator<T>> sketch2vec() const {
        std::vector<T, common::Allocator<T>> ret;
        auto &s = sketch();
        ret.reserve(s.size());
        for(const auto el: this->sketch()) ret.emplace_back(el);
        return ret;
    }
public:
    uint64_t sketch_size() const {
        return ss_;
    }
};




template<typename T, typename Cmp=std::greater<T>, typename Hasher=common::WangHash>
class KMinHash: public AbstractMinHash<T, Cmp> {
    std::vector<uint64_t, common::Allocator<uint64_t>> seeds_;
    std::vector<T, common::Allocator<uint64_t>> hashes_;
    Hasher hf_;
    // Uses k hash functions with seeds.
    // TODO: reuse code from hll/bf for vectorized hashing.
public:
    KMinHash(size_t nkeys, size_t sketch_size, uint64_t seedseed=137, Hasher &&hf=Hasher()):
        AbstractMinHash<T, Cmp>(sketch_size),
        hashes_(nkeys),
        hf_(std::move(hf))
    {
        DefaultRNGType gen(seedseed);
        seeds_.reserve(nseeds());
        while(seeds_.size() < nseeds()) seeds_.emplace_back(gen());
        throw NotImplementedError();
    }
    size_t nseeds() const {return this->nvals() * sizeof(uint64_t) / sizeof(T);}
    size_t nhashes_per_seed() const {return sizeof(uint64_t) / sizeof(T);}
    size_t nhashes_per_vector() const {return sizeof(uint64_t) / sizeof(Space::Type);}
    void addh(T val) {
        throw NotImplementedError();
        // Hash stuff.
        // Then use a vectorized minimization comparison instruction
    }
    SET_SKETCH(hashes_)
};

template<typename T, typename Allocator> struct FinalRMinHash; // Forward definition
/*
The sketch is the set of minimizers.

*/
template<typename T,
         typename Cmp=std::greater<T>,
         typename Hasher=common::WangHash,
         typename Allocator=common::Allocator<T>
        >
struct RangeMinHash: public AbstractMinHash<T, Cmp> {
protected:
    Hasher hf_;
    Cmp cmp_;
    ///using HeapType = std::priority_queue<T, std::vector<T, Allocator<T>>>;
    std::set<T, Cmp> minimizers_; // using std::greater<T> so that we can erase from begin()
    int counter; // Tracks the current size (number of elements) in the HashMap. Should be replaced with the size function when available.

public:
    using final_type = FinalRMinHash<T, Allocator>;
    using Compare = Cmp;
    RangeMinHash(size_t sketch_size, Hasher &&hf=Hasher(), Cmp &&cmp=Cmp()):
        AbstractMinHash<T, Cmp>(sketch_size), hf_(std::move(hf)), cmp_(std::move(cmp))
    {
    }
    ~RangeMinHash() {
    }
    RangeMinHash(std::string) {throw NotImplementedError("");}    
    /*
     * Using BCL
     * gets the maximum element from the set
     */
    T max_element() const {
        // fprintf(stderr, "RangeMinHash.max_element()\n");
#if 0
        for(const auto e: *this)
            assert(*begin() >= e);
#endif
        return *begin(); // NOTE: when begin() is updated, max_element() should be too. Just confirm that it's sorted properly
    }
    INLINE void addh(T val) {
        val = hf_(val);
        this->add(val);
    }
    INLINE void add(T val) {        
        if(minimizers_.size() == this->ss_) {
            if(cmp_(max_element(), val)) {
                minimizers_.insert(val);
                if(minimizers_.size() > this->ss_) {
                    minimizers_.erase(minimizers_.begin());
                }
            }
        } else minimizers_.insert(val);
    }
    auto begin() { fprintf(stderr, "RangeMinHash.begin()\n"); return minimizers_.begin();}
    auto begin() const {
        return minimizers_.begin();
    }
    auto end() { fprintf(stderr, "RangeMinHash.end()\n"); return minimizers_.end();}
    auto end() const {
        return minimizers_.end();
    }
    template<typename C2>
    size_t intersection_size(const C2 &o) const {
        return common::intersection_size(o, *this, cmp_);
    }
    double jaccard_index(const RangeMinHash &o) const {
        assert(minimizers_.size() == size() && size() == this->counter); // TODO: test. Just need to replace minimizers_.size() with new size
        double is = this->intersection_size(o);
        return is / (minimizers_.size() + o.size() - is);
    }
    final_type cfinalize() const {
        // fprintf(stderr, "RangeMinHash.cfinalize()\n");
        // std::vector<T> reta(minimizers_.begin(), minimizers_.end());
        // if(reta.size() < this->ss_)
        //     reta.insert(reta.end(), this->ss_ - reta.size(), std::numeric_limits<uint64_t>::max());
        // return final_type(std::move(reta));

        // TODO: test new version -- needs to be fixed
        std::vector<T> reta(minimizers_.begin(), minimizers_.end());
        if(reta.size() < this->ss_)
            reta.insert(reta.end(), this->ss_ - reta.size(), std::numeric_limits<uint64_t>::max());
        return final_type(std::move(reta));
    }
    final_type finalize() & {
        // TODO: add to tests
        fprintf(stderr, "RangeMinHash.finalize()\n");
        return cfinalize();
    }
    final_type finalize() const & {
        // TODO: add to tests
        fprintf(stderr, "RangeMinHash.finalize() const &\n");
        return cfinalize();
    }
    final_type finalize() && {
        // TODO: add to tests
        fprintf(stderr, "RangeMinHash.finalize() &&\n");
        return static_cast<const RangeMinHash &&>(*this).finalize();
    }
    final_type finalize() const && {
        // TODO: add to tests
        fprintf(stderr, "RangeMinHash.finalize() const && \n");
        std::vector<T> reta(minimizers_.begin(), minimizers_.end());
        if(reta.size() < this->ss_) {
            reta.insert(reta.end(), this->ss_ - reta.size(), std::numeric_limits<uint64_t>::max());
        }
        final_type ret(std::move(reta));
        return ret;
    }
    size_t size() const {
        return minimizers_.size();
    }
    using key_compare = typename decltype(minimizers_)::key_compare;
    SET_SKETCH(minimizers_)
};

namespace weight {


struct EqualWeight {
    // Do you even weight bro?
    template<typename T>
    constexpr double operator()(T &x) const {return 1.;}
};

template<typename ScorerType>
struct DictWeight {
    const ScorerType &m_;
    DictWeight(const ScorerType &m): m_(m) {}
    template<typename T>
    double operator()(T &x) const {
        return m_.score(x);
    }
};
}


template<typename T>
double jaccard_index(const T &a, const T &b) {
    return a.jaccard_index(b);
}

} // inline namespace minhash
namespace mpi_mh = mpi_minhash;
} // namespace sketch
