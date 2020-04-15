#pragma once
#include <mutex>
//#include <queue>
#include "hll.h" // For common.h and clz functions
#include "fixed_vector.h"
#include <unordered_map>
#include "isz.h"
#include <bcl/bcl.hpp>
#include <bcl/containers/HashMap.hpp>
#include <bcl/containers/HashSet.hpp>


/*
 * TODO: support minhash using sketch size and a variable number of hashes.
 * Implementations: KMinHash (multiple hash functions)
 *                  RangeMinHash (the lowest range in a set, but only uses one hash)
 *                  HyperMinHash
 */

namespace sketch {
inline namespace bcl_minhash {

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
    BCL::HashMap<T, int> * bcl_map_; // Use of hash map is a temporary replacement for a set. When Sets are available in BCL.
    BCL::HashSet<T,int> * bcl_set_;

public:
    using final_type = FinalRMinHash<T, Allocator>;
    using Compare = Cmp;
    RangeMinHash(size_t sketch_size, Hasher &&hf=Hasher(), Cmp &&cmp=Cmp()):
        AbstractMinHash<T, Cmp>(sketch_size), hf_(std::move(hf)), cmp_(std::move(cmp))
    {
        bcl_set_ = new BCL::HashSet<T,int>(sketch_size);
        bcl_map_ = new BCL::HashMap<T,int>(sketch_size);
    }
    ~RangeMinHash() {
        delete bcl_map_;
        delete bcl_set_;
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
        // fprintf(stderr, "RangeMinHash.addh()\n");
        val = hf_(val);
        this->add(val);

        // NOTE: updating add() should update this as well
    }
    INLINE void add(T val) {
        // fprintf(stderr, "RangeMinHash.add()\n");
        
        // TODO: make work for BCL map
        // if(this->bcl_set_->size() == this->ss_) {
        //     if(cmp_(max_element(), val)) {
        //         bcl_set_->insert_or_assign(val, 1);
        //         if(bcl_set_->size() > this->ss_) { 
        //             bcl_set_->delete_max();
        //             // TODO: figure out what's going on
        //             // bcl_set_->insert_or_assign(*minimizers_.begin(), 0);
        //             // auto test = *bcl_map_->begin();
        //             // auto test2 = test.value_type();
        //             // auto test3 = test2.first;
        //             // bcl_map_->insert_or_assign((*bcl_map_->begin()), 0); // TODO: get info more cleanly
        //         }
        //     }
        // } else {
        //     bcl_set_->insert_or_assign(val, 1);
        // }

        // NOTE: test new version (above)
        // if(minimizers_.size() != bcl_map_->capacity()) {
        //     fprintf(stderr, "Testing size: min = %ld, bcl = %ld, bcl local = %ld\n", minimizers_.size(), bcl_map_->capacity(), bcl_map_->local_capacity());
        // }

        if(minimizers_.size() == this->ss_) {
            if(cmp_(max_element(), val)) {
                minimizers_.insert(val);
                bcl_set_->insert_or_assign(val, 1);
                if(minimizers_.size() > this->ss_) {
                    minimizers_.erase(minimizers_.begin());
                    // bcl_set_->delete_max(); //TODO: including delete_max() here causes an error, do not use in current state
                    // bcl_map_->insert_or_assign(minimizers_.begin()->val, 0);
                }
            }
        } else minimizers_.insert(val);
    }
    // TODO: add to tests and update
    auto begin() { fprintf(stderr, "RangeMinHash.begin()\n"); return minimizers_.begin();}
    // TODO: add to tests and update
    auto begin() const {
        // assert(minimizers_->begin() == bcl_map_->begin());
        // return bcl_map_->begin(); // TODO: use in future
        return minimizers_.begin();
    }
    // TODO: add to tests and update
    auto end() { fprintf(stderr, "RangeMinHash.end()\n"); return minimizers_.end();}
    auto end() const {
        // assert(minimizers_->end() == bcl_map_->end());
        // return bcl_map_->end(); // TODO: use in future
        return minimizers_.end();
    }
    template<typename C2>
    size_t intersection_size(const C2 &o) const {
        // fprintf(stderr, "RangeMinHash.intersection_size()\n");
        return common::intersection_size(o, *this, cmp_);

        // TODO: figure out way to use intersection size with this new version???
    }
    double jaccard_index(const RangeMinHash &o) const {
        assert(minimizers_.size() == size() && size() == this->bcl_set_->size()); // TODO: test. Just need to replace minimizers_.size() with new size
        // fprintf(stderr, "RangeMinHash.jaccard_index()\n");
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
        // std::vector<T> reta(bcl_map_->begin(), bcl_map_->end()); // TODO: make this line work. Right now, a vector can't be built from a BCL::GlobalhashMapIterator
        if(reta.size() < this->ss_)
            reta.insert(reta.end(), this->ss_ - reta.size(), std::numeric_limits<uint64_t>::max());
        return final_type(std::move(reta));
        // TODO: does the vector need to be a BCL structure?

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
        if(minimizers_.size() != bcl_set_->size()) { // TODO: delete debug print statement
            fprintf(stderr, "Testing size: min = %ld, bcl = %ld, size = %ld\n", minimizers_.size(), bcl_map_->capacity(), this->bcl_set_->size());
        }
        // if(minimizers_.size() != bcl_map_->capacity()) {
        //     fprintf(stderr, "Testing size: min = %ld, bcl = %ld, bcl local = %ld\n", minimizers_.size(), bcl_map_->capacity(), bcl_map_->local_capacity());
        // }
        return bcl_set_->size();
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

template<typename T, typename Allocator=std::allocator<T>>
struct FinalRMinHash {
    static_assert(std::is_unsigned<T>::value, "must be unsigned btw");

    using allocator = Allocator;
    std::vector<T, Allocator> first;
    uint64_t lsum_;
    using container_type = decltype(first);
    size_t intersection_size(const FinalRMinHash &o) const {
        return common::intersection_size(first, o.first);
    }
    double jaccard_index(const FinalRMinHash &o) const {
        double is = intersection_size(o);
        return is / ((size() << 1) - is);
    }
    void print_all_cards() const {
        for(const auto x: {HARMONIC_MEAN, ARITHMETIC_MEAN, MEDIAN}) {
            std::fprintf(stderr, "cardest with x = %d is %lf\n", int(x), cardinality_estimate(x));
        }
    }
    FinalRMinHash &operator+=(const FinalRMinHash &o) {
        std::vector<T, Allocator> newfirst; newfirst.reserve(o.size());
        if(this->size() != o.size()) throw std::runtime_error("Non-matching parameters for FinalRMinHash comparison");
        auto i1 = this->begin(), i2 = o.begin();
        while(newfirst.size() < first.size()) {
            if(*i2 < *i1) {
                newfirst.push_back(*i2++);
            } else if(*i1 < *i2) {
                newfirst.push_back(*i1++);
            } else newfirst.push_back(*i1), ++i1, ++i2;
        }
        std::swap(newfirst, first);
        return *this;
    }
    FinalRMinHash operator+(const FinalRMinHash &o) const {
        auto tmp = *this;
        tmp += o;
        return tmp;
    }
    double union_size(const FinalRMinHash &o) const {
        PREC_REQ(this->size() == o.size(), "Non-matching parameters for FinalRMinHash comparison");
        size_t n_in_sketch = 0;
        auto i1 = this->rbegin(), i2 = o.rbegin();
        T mv;
        while(n_in_sketch < first.size() - 1) {
            // Easier to branch-predict:  http://www.vldb.org/pvldb/vol8/p293-inoue.pdf
            if(*i1 != *i2) ++i1, ++i2;
            else {
                const int c = *i1 < *i2;
                i2 += !c; i1 += c;
            }
            ++n_in_sketch;
        }
        mv = *i1 < *i2 ? *i1: *i2;
        // TODO: test after refactoring
        assert(i1 < this->rend());
        return double(std::numeric_limits<T>::max()) / (mv) * this->size();
    }
    double cardinality_estimate(MHCardinalityMode mode=ARITHMETIC_MEAN) const {
        // KMV estimate
        double sum = (std::numeric_limits<T>::max() / double(this->max_element()) * first.size());
        return sum;
    }
    void sum() const {
        lsum_ = std::accumulate(this->first.begin(), this->first.end(), size_t(0));
    }
    template<typename WeightFn=weight::EqualWeight>
    double tf_idf(const FinalRMinHash &o, const WeightFn &fn=WeightFn()) const {
        const size_t lsz = size();
        const size_t rsz = o.size();
        double num = 0, denom = 0;
        size_t i1 = 0, i2 = 0;
        for(;;) {
            if(first[i1] < o.first[i2]) {
                denom += fn(first[i1]);
                if(++i1 == lsz) break;
            } else if(o.first[i2] < first[i1]) {
                denom += fn(o.first[i2]);
                assert(i2 != rsz);
                if(++i2 == rsz) break;
            } else {
                const auto v1 = fn(first[i1]), v2 = fn(o.first[i2]);
                denom += std::max(v1, v2);
                num += std::min(v1, v2);
                assert(i2 != rsz);
                ++i1; ++i2;
                if(i1 == lsz || i2 == rsz) break;
            }
        }
        while(i1 < lsz) denom += fn(first[i1++]);
        while(i2 < rsz) denom += fn(o.first[i2++]);
        return num / denom;
    }
    typename container_type::const_iterator begin() {
        return first.cbegin();
    }
    typename container_type::const_iterator end() {
        return first.cend();
    }
    typename container_type::const_iterator begin() const {
        return first.cbegin();
    }
    typename container_type::const_iterator end() const {
        return first.cend();
    }
    typename container_type::const_reverse_iterator rbegin() {
        return first.crbegin();
    }
    typename container_type::const_reverse_iterator rend() {
        return first.crend();
    }
    typename container_type::const_reverse_iterator rbegin() const {
        return first.crbegin();
    }
    typename container_type::const_reverse_iterator rend() const {
        return first.crend();
    }
    void free() {
        decltype(first) tmp; std::swap(tmp, first);
    }
    template<typename Alloc>
    FinalRMinHash(const std::vector<T, Alloc> &ofirst): first(ofirst.size()) {std::copy(ofirst.begin(), ofirst.end(), first.begin()); sort();}
    template<typename It>
    FinalRMinHash(It start, It end): first(std::distance(start, end)) {
        std::copy(start, end, first.begin());
        sort();
    }
    template<typename Alloc, typename=std::enable_if_t<std::is_same<Alloc, allocator>::value>>
    FinalRMinHash(std::vector<T, Alloc> &&ofirst): first(std::move(ofirst)) {
        sort();
    }
    FinalRMinHash(FinalRMinHash &&o): first(std::move(o.first)) {sort();}
    ssize_t read(gzFile fp) {
        uint64_t sz;
        if(gzread(fp, &sz, sizeof(sz)) != sizeof(sz)) throw ZlibError("Failed to read");
        first.resize(sz);
        ssize_t ret = sizeof(sz), nb = first.size() * sizeof(first[0]);
        if(gzread(fp, first.data(), nb) != nb) throw ZlibError("Failed to read");
        ret += nb;
        sort();
        return ret;
    }
    ssize_t write(gzFile fp) const {
        uint64_t sz = this->size();
        ssize_t ret = gzwrite(fp, &sz, sizeof(sz));
        ret += gzwrite(fp, first.data(), first.size() * sizeof(first[0]));
        return ret;
    }
    DBSKETCH_READ_STRING_MACROS
    DBSKETCH_WRITE_STRING_MACROS
    auto max_element() const {
        const auto ret = first.back();
        assert(std::accumulate(first.begin(), first.end(), true, [&](bool t, auto v) {return t && ret >= v;}));
        return ret;
    }
    template<typename Hasher, typename Cmp>
    FinalRMinHash(RangeMinHash<T, Cmp, Hasher> &&prefinal): FinalRMinHash(std::move(prefinal.finalize())) {
        prefinal.clear();
        sort();
    }
    FinalRMinHash(const std::string &s): FinalRMinHash(s.data()) {}
    FinalRMinHash(gzFile fp) {read(fp);}
    FinalRMinHash(const char *infname) {read(infname);}
    size_t size() const {return first.size();}
protected:
    FinalRMinHash() {}
    FinalRMinHash(const FinalRMinHash &o) = default;
    FinalRMinHash &operator=(const FinalRMinHash &o) = default;
    void sort() {
        common::sort::default_sort(this->first.begin(), this->first.end());
    }
};

template<typename T, typename CountType> struct FinalCRMinHash; // Forward


template<typename T,
         typename Cmp=std::greater<T>,
         typename Hasher=common::WangHash,
         typename CountType=uint32_t
        >
class CountingRangeMinHash: public AbstractMinHash<T, Cmp> {
    static_assert(std::is_arithmetic<CountType>::value, "CountType must be arithmetic");
    struct VType {
        T first;
        mutable CountType second;
        inline bool operator<(const VType &b) const {
            return Cmp()(this->first , b.first);
        }
        inline bool operator>(const VType &b) const {
            return !Cmp()(this->first , b.first);
        }
        inline bool operator==(const VType &b) const {
            return this->first == b.first;
        }
        VType(T v, CountType c): first(v), second(c) {}
        VType(const VType &o): first(o.first), second(o.second) {}
        VType(gzFile fp) {if(gzread(fp, this, sizeof(*this)) != sizeof(*this)) throw ZlibError("Failed to read");}
    };
    Hasher hf_;
    Cmp cmp_;
    mutable CountType cached_sum_sq_ = 0, cached_sum_ = 0;
    std::set<VType> minimizers_; // using std::greater<T> so that we can erase from begin()
public:
    const auto &min() const {return minimizers_;}
    using size_type = CountType;
    using final_type = FinalCRMinHash<T, CountType>;
    auto size() const {return minimizers_.size();}
    auto begin() {return minimizers_.begin();}
    auto begin() const {return minimizers_.begin();}
    auto rbegin() {return minimizers_.rbegin();}
    auto rbegin() const {return minimizers_.rbegin();}
    auto end() {return minimizers_.end();}
    auto end() const {return minimizers_.end();}
    auto rend() {return minimizers_.rend();}
    auto rend() const {return minimizers_.rend();}
    void free() {
        std::set<VType> tmp;
        std::swap(tmp, minimizers_);
    }
    CountingRangeMinHash(size_t n, Hasher &&hf=Hasher(), Cmp &&cmp=Cmp()): AbstractMinHash<T, Cmp>(n), hf_(std::move(hf)), cmp_(std::move(cmp)) {}
    CountingRangeMinHash(std::string s): CountingRangeMinHash(0) {throw NotImplementedError("");}
    double cardinality_estimate(MHCardinalityMode mode=ARITHMETIC_MEAN) const {
        return double(std::numeric_limits<T>::max()) / std::max_element(minimizers_.begin(), minimizers_.end(), [](auto x, auto y) {return x.first < y.first;})->first * minimizers_.size();
    }
    INLINE void add(T val) {
        if(minimizers_.size() == this->ss_) {
            if(cmp_(begin()->first, val)) {
                auto it = minimizers_.find(VType(val, 0));
                if(it == minimizers_.end()) {
                    minimizers_.erase(begin());
                    minimizers_.insert(VType(val, CountType(1)));
                } else ++it->second;
            }
        } else minimizers_.insert(VType(val, CountType(1)));
    }
    INLINE void addh(T val) {
        val = hf_(val);
        this->add(val);
    }
    auto max_element() const {
        return minimizers_.begin()->first;
    }
    auto sum_sq() {
        if(cached_sum_sq_) goto end;
        cached_sum_sq_ = std::accumulate(this->begin(), this->end(), 0., [](auto s, const VType &x) {return s + x.second * x.second;});
        end:
        return cached_sum_sq_;
    }
    auto sum_sq() const {return cached_sum_sq_;}
    auto sum() {
        if(cached_sum_) goto end;
        cached_sum_ = std::accumulate(std::next(this->begin()), this->end(), this->begin()->second, [](auto s, const VType &x) {return s + x.second;});
        end:
        return cached_sum_;
    }
    auto sum() const {return cached_sum_;}
    double histogram_intersection(const CountingRangeMinHash &o) const {
        assert(o.size() == size());
        size_t denom = 0, num = 0;
        auto i1 = minimizers_.begin(), i2 = o.minimizers_.begin();
        auto e1 = minimizers_.end(), e2 = o.minimizers_.end();
        for(;;) {
            if(cmp_(i1->first, i2->first)) {
                denom += (i1++)->second;
                if(i1 == e1) break;
            } else if(cmp_(i2->first, i1->first)) {
                denom += (i2++)->second;
                if(i2 == e2) break;
            } else {
                const auto v1 = i1->second, v2 = i2->second;
                denom += std::max(v1, v2);
                num += std::min(v1, v2);
                ++i1; ++i2;
                if(i1 == e1) break;
                if(i2 == e2) break;
            }
        }
        while(i1 != e1) denom += i1++->second;
        while(i2 != e2) denom += i2++->second;
        return static_cast<double>(num) / denom;
    }
    double containment_index(const CountingRangeMinHash &o) const {
        assert(o.size() == size());
        size_t denom = 0, num = 0;
        auto i1 = minimizers_.begin(), i2 = o.minimizers_.begin(),
             e1 = minimizers_.end(),   e2 = o.minimizers_.end();
        for(;;) {
            if(cmp_(i1->first, i2->first)) {
                denom += i1->second;
                if(++i1 == e1) break;
            } else if(cmp_(i2->first, i1->first)) {
                if(++i2 == e2) break;
            } else {
                const auto v1 = i1->second, v2 = i2->second;
                denom += v1;
                num += std::min(v1, v2);
                if(++i1 == e1) break;
                if(++i2 == e2) break;
            }
        }
        return static_cast<double>(num) / denom;
    }
    DBSKETCH_WRITE_STRING_MACROS
    DBSKETCH_READ_STRING_MACROS
    ssize_t write(gzFile fp) const {
        uint64_t n = minimizers_.size();
        if(gzwrite(fp, &n, sizeof(n)) != sizeof(n)) throw ZlibError("Failed to write");
        for(const auto &pair: minimizers_) {
            if(gzwrite(fp, std::addressof(pair), sizeof(pair)) != sizeof(pair))
                ZlibError("Failed to write");
        }
        return sizeof(VType) * minimizers_.size();
    }
    ssize_t read(gzFile fp) {
        uint64_t n;
        if(gzread(fp, &n, sizeof(n)) != sizeof(n)) throw ZlibError("Failed to read");
        for(size_t i = n; i--; minimizers_.insert(VType(fp)));
        return sizeof(n) + sizeof(VType) * n;
    }

    void clear() {
        decltype(minimizers_) tmp;
        std::swap(tmp, minimizers_);
    }
    template<typename WeightFn=weight::EqualWeight>
    double tf_idf(const CountingRangeMinHash &o, const WeightFn &fn) const {
        assert(o.size() == size());
        double denom = 0, num = 0;
        auto i1 = minimizers_.begin(), i2 = o.minimizers_.begin();
        const auto e1 = minimizers_.cend(), e2 = o.minimizers_.cend();
        for(;;) {
            auto &lhf = i1->first, &rhf = i2->first;
            if(cmp_(lhf, rhf)) {
                denom += (i1->second) * fn(lhf);
                if(++i1 == e1) break;
            } else if(cmp_(rhf, lhf)) {
                denom += (i2->second) * fn(rhf);
                if(++i2 == e2) break;
            } else {
                assert(rhf == lhf);
                const auto v1 = i1->second * fn(lhf), v2 = i2->second * fn(rhf);
                denom += std::max(v1, v2);
                num += std::min(v1, v2);
                if(++i1 == e1 || ++i2 == e2) break;
            }
        }
        while(i1 < e1) denom += i1->second * fn(i1->first), ++i1;
        while(i2 < e2) denom += i2->second * fn(i2->first), ++i2;
        return static_cast<double>(num) / denom;
    }
    final_type finalize() & {
        return cfinalize();
    }
    final_type finalize() const & {
        return cfinalize();
    }
    final_type cfinalize() const & {
        return FinalCRMinHash<T, CountType>(*this);
    }
    final_type finalize() && {
        return static_cast<const CountingRangeMinHash &&>(*this).finalize();
    }
    final_type finalize() const && {
        auto ret(FinalCRMinHash<T, CountType>(std::move(*this)));
        return ret;
    }
    template<typename Func>
    void for_each(const Func &func) const {
        for(const auto &i: minimizers_) {
            func(i);
        }
    }
    template<typename Func>
    void for_each(const Func &func) {
        for(auto &i: minimizers_) {
            func(i);
        }
    }
    void print() const {
        for_each([](auto &p) {std::fprintf(stderr, "key %s with value %zu\n", std::to_string(p.first).data(), size_t(p.second));});
    }
    size_t union_size(const CountingRangeMinHash &o) const {
        size_t denom = 0;
        auto i1 = minimizers_.begin(), i2 = o.minimizers_.begin();
        auto e1 = minimizers_.end(), e2 = o.minimizers_.end();
        for(;;) {
            if(cmp_(i1->first, i2->first)) {
                denom += i1->second;
                if(++i1 == e1) break;
            } else if(cmp_(i2->first, i1->first)) {
                denom += i2->second;
                if(++i2 == e2) break;
            } else {
                denom += std::max(i1->second, i2->second);
                if( (++i1 == e1) | (++i2 == e2)) break;
            }
        }
        while(i1 != e1) denom += i1->second, ++i1;
        while(i2 != e2) denom += i2->second, ++i2;
        return denom;
    }
    size_t intersection_size(const CountingRangeMinHash &o) const {
        size_t num = 0;
        auto i1 = minimizers_.begin(), i2 = o.minimizers_.begin();
        auto e1 = minimizers_.end(), e2 = o.minimizers_.end();
        for(;;) {
            if(cmp_(i1->first, i2->first)) {
                if(++i1 == e1) break;
            } else if(cmp_(i2->first, i1->first)) {
                if(++i2 == e2) break;
            } else {
                num += std::min(i1->second, i2->second);
                if(++i1 == e1) break;
                if(++i2 == e2) break;
            }
        }
        return num;
    }
    template<typename C2>
    double jaccard_index(const C2 &o) const {
        double is = this->intersection_size(o);
        return is / (minimizers_.size() + o.size() - is);
    }
};



template<typename T, typename CountType=uint32_t>
struct FinalCRMinHash: public FinalRMinHash<T> {
    using super = FinalRMinHash<T>;
    std::vector<CountType> second;
    uint64_t count_sum_;
    double count_sum_l2norm_;
    using count_type = CountType;
    using key_type = T;
    FinalCRMinHash(gzFile fp) {
        this->read(fp);
    }
    FinalCRMinHash(const std::string &path): FinalCRMinHash(path.data()) {}
    FinalCRMinHash(const char *path) {
        this->read(path);
    }
    void free() {
        super::free();
        std::vector<CountType> tmp;
        std::swap(tmp, second);
    }
    size_t countsum() const {return std::accumulate(second.begin(), second.end(), size_t(0), [](auto sz, auto sz2) {return sz += sz2;});}
    size_t countsumsq() const {
#if !NDEBUG
        size_t lsum = 0;
        for(const auto v: this->second) lsum += v * v;
        
        size_t asum = std::accumulate(second.begin(), second.end(), size_t(0), [](auto sz, auto sz2) {return sz += sz2 * sz2;});
        assert(asum == lsum);
        return asum;
#else
        return std::accumulate(second.begin(), second.end(), size_t(0), [](auto sz, auto sz2) {return sz += sz2 * sz2;});
#endif
    }
    double cosine_distance(const FinalCRMinHash &o) const {
        return dot(o) / count_sum_l2norm_ / o.count_sum_l2norm_;
    }
    double dot(const FinalCRMinHash &o) const {
        const size_t lsz = this->size(), rsz = o.size();
        size_t num = 0;
        for(size_t i1 = 0, i2 = 0;;) {
            if(this->first[i1] < o.first[i2]) {
                if(++i1 == lsz) break;
            } else if(o.first[i2] < this->first[i1]) {
                if(++i2 == rsz) break;
            } else {
                const auto v1 = o.second[i2], v2 = second[i1];
                num += v1 * v2;
                if(++i1 == lsz || ++i2 == rsz) break;
            }
        }
        return static_cast<double>(num);
    }
    double histogram_intersection(const FinalCRMinHash &o) const {
        const size_t lsz = this->size(), rsz = o.size();
        assert(std::accumulate(this->second.begin(), this->second.end(), size_t(0)) == this->count_sum_);
        assert(std::accumulate(o.second.begin(), o.second.end(), size_t(0)) == o.count_sum_);
        assert(std::sqrt(std::accumulate(this->second.begin(), this->second.end(), 0., [](auto psum, auto newv) {return psum + newv * newv;})) == this->count_sum_l2norm_);
        assert(std::sqrt(std::accumulate(o.second.begin(), o.second.end(), 0., [](auto psum, auto newv) {return psum + newv * newv;})) == o.count_sum_l2norm_);
        assert(count_sum_ > 0);
        assert(o.count_sum_ > 0);
        size_t num = 0;
        size_t i1 = 0, i2 = 0;
        for(;;) {
            auto &lhs = this->first[i1];
            auto &rhs = o.first[i2];
            if(lhs < rhs) {
                if(++i1 == lsz) break;
            } else if(rhs < lhs) {
                if(++i2 == rsz) break;
            } else {
                assert(!(lhs < rhs));
                assert(!(rhs < lhs));
                assert(lhs == rhs);
                auto lhv = this->second[i1++], rhv = o.second[i2++];
                num += std::min(lhv, rhv);
                if(i1 == lsz || i2 == rsz) break;
            }
        }
        return static_cast<double>(num) / (count_sum_ + o.count_sum_ - num);
    }
    DBSKETCH_READ_STRING_MACROS
    DBSKETCH_WRITE_STRING_MACROS
    ssize_t read(gzFile fp) {
        uint64_t nelem;
        ssize_t ret = gzread(fp, &nelem, sizeof(nelem));
        this->first.resize(nelem);
        ret += gzread(fp, &count_sum_, sizeof(count_sum_));
        ret += gzread(fp, &count_sum_l2norm_, sizeof(count_sum_l2norm_));
        ret += gzread(fp, this->first.data(), sizeof(this->first[0]) * nelem);
        this->second.resize(nelem);
        ret += gzread(fp, this->second.data(), sizeof(this->second[0]) * nelem);
        assert(std::set<uint64_t>(this->first.begin(), this->first.end()).size() == this->first.size());
        prepare();
        return ret;
    }
    double union_size(const FinalCRMinHash &o) const {
        PREC_REQ(this->size() == o.size(), "mismatched parameters");
        return super::union_size(o);
    }
    ssize_t write(gzFile fp) const {
        assert(std::set<uint64_t>(this->first.begin(), this->first.end()).size() == this->first.size());
        uint64_t nelem = second.size();
        ssize_t ret = gzwrite(fp, &nelem, sizeof(nelem));
        ret += gzwrite(fp, &count_sum_, sizeof(count_sum_));
        ret += gzwrite(fp, &count_sum_l2norm_, sizeof(count_sum_l2norm_));
        ret += gzwrite(fp, this->first.data(), sizeof(this->first[0]) * nelem);
        ret += gzwrite(fp, this->second.data(), sizeof(this->second[0]) * nelem);
        return ret;
    }
    template<typename WeightFn=weight::EqualWeight>
    double tf_idf(const FinalCRMinHash &o, const WeightFn &fn=WeightFn()) const {
        const size_t lsz = this->size();
        const size_t rsz = o.size();
        assert(rsz == o.second.size());
        assert(rsz == o.first.size());
        double denom = 0, num = 0;
        size_t i1 = 0, i2 = 0;
        for(;;) {
            auto &lhs = this->first[i1];
            auto &rhs = o.first[i2];
            const auto lhv = second[i1] * fn(lhs), rhv = o.second[i2] * fn(rhs);
            if(lhs < rhs) {
                denom += lhv * fn(lhs);
                if(++i1 == lsz) break;
            } else if(rhs < lhs) {
                denom += rhv * fn(rhs);
                if(++i2 == rsz) break;
            } else {
                denom += std::max(lhv, rhv);
                num += std::min(lhv, rhv);
                ++i2, ++i1;
                if(i2 == rsz || i1 == lsz) break;
            }
            assert(i2 < o.second.size());
        }
        while(i1 < lsz) {
            denom += this->second.operator[](i1) * fn(this->first.operator[](i1)), ++i1;
            //denom += this->first[i1] * fn(second[i1]), ++i1;
        }
        while(i2 < rsz) {
            denom += o.second.operator[](i2) * fn(o.first.operator[](i2)), ++i2;
        }
        //std::fprintf(stderr, "[%s] after finishing num %f, denom %f\n", __PRETTY_FUNCTION__, num, denom);
        return num / denom;
    }
    void prepare(size_t ss=0) {
        const size_t fs = this->first.size();
        ss = ss ? ss: fs;
        fixed::vector<std::pair<key_type, count_type>> tmp(fs);
        for(size_t i = 0; i < fs; ++i)
            tmp[i] = std::make_pair(this->first[i], second[i]);
        auto tmpcmp = [](auto x, auto y) {return x.first < y.first;};
        common::sort::default_sort(tmp.begin(), tmp.end(), tmpcmp);
        for(size_t i = 0; i < this->first.size(); ++i)
            this->first[i] = tmp[i].first, this->second[i] = tmp[i].second;
        assert(std::is_sorted(this->first.begin(), this->first.end()));
        PREC_REQ(this->first.size() == this->second.size(), "Counts and hashes must have equal size");
        const std::ptrdiff_t diff = ss - this->first.size();
        if(diff > 0) {
            this->first.insert(this->first.end(), diff, std::numeric_limits<T>::max());
            this->second.insert(this->second.end(), diff, 0);
        } else if(diff < 0) {
            const auto fe = this->first.end();
            const auto se = this->second.end();
            assert(fe + diff < fe);
            assert(se + diff < se);
            this->first.erase(fe + diff, fe);
            this->second.erase(se + diff, se);
        }
        assert(this->first.size() == this->second.size());
        assert(this->first.size() == ss);
        count_sum_ = countsum();
        count_sum_l2norm_ = std::sqrt(countsumsq());
        POST_REQ(this->first.size() == this->second.size(), "Counts and hashes must have equal size");
    }
    template<typename Valloc, typename=std::enable_if_t<!std::is_same<Valloc, typename super::allocator>::value>>
    FinalCRMinHash(std::vector<T, Valloc> &&first, std::vector<CountType> &&second, size_t ss=0) {
        //this->first = std::move(first);
        this->first.assign(first.begin(), first.end());
        this->second = std::move(second);
        std::swap(first, std::vector<T, Valloc>());
        prepare(ss);
    }
    FinalCRMinHash(std::vector<T, typename super::allocator> &&first, std::vector<CountType> &&second, size_t ss=0) {
#if VERBOSE_AF
        std::fprintf(stderr, "first size: %zu.\n", first.size());
        std::fprintf(stderr, "second size: %zu.\n", second.size());
#endif
        this->first = std::move(first);
        this->second = std::move(second);
        prepare(ss);
    }
    template<typename Valloc>
    FinalCRMinHash(std::pair<std::vector<T, Valloc>, std::vector<CountType>> &&args, size_t ss=0):
        FinalCRMinHash(std::move(args.first), std::move(args.second), ss)
    {
    }
    template<typename Hasher, typename Cmp>
    static std::pair<std::vector<T>, std::vector<CountType>> crmh2vecs(const CountingRangeMinHash<T, Cmp, Hasher, CountType> &prefinal) {
        std::vector<T> tmp;
        std::vector<CountType> tmp2;
        tmp.reserve(prefinal.size()), tmp2.reserve(prefinal.size());
        for(const auto &pair: prefinal)
            tmp.push_back(pair.first), tmp2.push_back(pair.second);
        return std::make_pair(std::move(tmp), std::move(tmp2));
    }
    template<typename Hasher, typename Cmp>
    FinalCRMinHash(const CountingRangeMinHash<T, Cmp, Hasher, CountType> &prefinal): FinalCRMinHash(crmh2vecs(prefinal), prefinal.sketch_size()) {
        assert(this->first.size() == prefinal.sketch_size());
    }
    template<typename Hasher, typename Cmp>
    FinalCRMinHash(CountingRangeMinHash<T, Cmp, Hasher, CountType> &&prefinal): FinalCRMinHash(static_cast<const CountingRangeMinHash<T, Cmp, Hasher, CountType> &>(prefinal)) {
        prefinal.clear();
    }
    double jaccard_index(const FinalCRMinHash &o) const {
        double us = union_size(o);
        double sz1 = cardinality_estimate(ARITHMETIC_MEAN);
        double sz2 = o.cardinality_estimate(ARITHMETIC_MEAN);
        double is = sz1 + sz2 - us;
        return std::max(0., is / us);
    }
    double containment_index(const FinalCRMinHash &o) const {
        double us = union_size(o);
        double sz1 = cardinality_estimate();
        double sz2 = o.cardinality_estimate();
        double is = sz1 + sz2 - us;
        return std::max(0., is / us);
    }
    double cardinality_estimate(MHCardinalityMode mode=ARITHMETIC_MEAN) const {
        return FinalRMinHash<T>::cardinality_estimate(mode);
    }
};

template<typename T>
double jaccard_index(const T &a, const T &b) {
    return a.jaccard_index(b);
}

} // inline namespace minhash
namespace bcl_mh = bcl_minhash;
} // namespace sketch
