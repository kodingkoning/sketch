#pragma once
#include "blaze/Math.h"
#include "div.h"
#include "common.h"
#include <queue>
#if 0
#if __cplusplus <= 201703L
#error("Need c++20")
#endif
#include <concepts>
#endif

namespace sketch {

inline namespace rnla {
enum Transform {
    CountSketch = 0,
    WoodruffZhang = 1
};

/***
   Utilities for compressing and decompressing
  */

//template<typename FT, bool SO>
//blaze::Dynamic

// TODO: write a generic method compatible with pytorch/ATen
//       (this actually does work, assuming we're working on the last level of a packedtensoraccessor), just not on GPU yet.
//       1. Random Laplace Transforms (Random Laplace Feature Maps for Semigroup Kernels on Histograms, Quasi-Monte Carlo Feature Maps for Shift-Invariant Kernels)
//       2. Tensor sketch: Count Sketch, but perform FFT on inputs and IFFT on outputs
//       3. Quasi Random Fourier Transform
//       4. Full matrix processing
//

struct CountApplicator {
    template<typename VT, typename RNG=blaze::RNG>
    VT apply(VT v, uint64_t rem) const {
        return v * (rem & 1 ? 1: -1);
    }
    template<typename VT, typename RNG=blaze::RNG>
    auto inverse(VT v, uint64_t rem) const {
        return v * (rem & 1 ? 1: -1);
    }
};

template<typename Distribution>
struct DistributionalApplicator {
    mutable Distribution dist_;

    using param_type = typename Distribution::param_type;
    template<typename...Args>
    DistributionalApplicator(Args &&...args): dist_(std::forward<Args>(args)...) {}

    template<typename VT, typename RNG=blaze::RNG>
    VT apply(VT val, uint64_t rem) {
        wy::WyHash<uint64_t, 1> rng(rem >> 1);
        auto mult = dist_(rng) * (rem & 1 ? 1.: -1.);
        return val * mult;
    }

    template<typename VT, typename RNG=blaze::RNG>
    auto inverse(VT v, uint64_t rem) {
        wy::WyHash<uint64_t, 1> rng(rem >> 1);
        auto mult = dist_(rng) * (rem & 1 ? 1.: -1.);
        return v / mult;
    }
};
struct WZApplicator: public DistributionalApplicator<std::exponential_distribution<double>> {};


template<typename C, typename C2, typename Hasher, typename Applicator>
auto &generic_compress(const C &in, C2 &ret, size_t newdim, const Hasher &hf, const Applicator &app) {
    if(newdim > in.size()) throw 1;
    std::fill(ret.begin(), ret.end(), static_cast<std::decay_t<decltype(ret[0])>>(0));
    PREC_REQ(newdim <= in.size(), "newdim cannot be larger");
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(newdim);
    PREC_REQ(ret.size() == newdim * ns, "out size doesn't match parameters");
    for(unsigned j = 0; j < ns; ++j) {
        for(unsigned i = 0; i < in.size(); ++i) {
            auto v = in[i];
            auto hv = hf(i, j);
            auto dm = div.divmod(hv);
            auto ind = dm.rem * ns + j;
            ret[ind] += app.apply(v, dm.quot);
        }
    }
    return ret;
}

template<typename FT, typename C2, bool SO, typename Hasher=KWiseHasherSet<4>, typename Applicator=CountApplicator>
auto &generic_compress(const blaze::CompressedVector<FT, SO> &in, C2 &ret, size_t newdim, const Hasher &hf=Hasher(), const Applicator &applicator=Applicator()) {
    if(newdim > in.size()) throw 1;
    std::fill(ret.begin(), ret.end(), static_cast<std::decay_t<decltype(ret[0])>>(0));
    PREC_REQ(newdim <= in.size(), "newdim cannot be larger");
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(newdim);
    PREC_REQ(ret.size() == newdim * ns, "out size doesn't match parameters");
    for(const auto &pair: in) {
        const auto idx = pair.index();
        const auto v = pair.value();
        for(unsigned j = 0; j < ns; ++j) {
            auto hv = hf(idx, j);
            auto dm = div.divmod(hv);
            auto ind = dm.rem * ns + j;
            ret[ind] += applicator.apply(v, dm.quot);
        }
    }
    return ret;
}

template<typename C, typename Hasher=KWiseHasherSet<4>, typename Applicator=CountApplicator>
auto generic_compress(const C &in, size_t newdim, const Hasher &hf, const Applicator &applicator) {
    C ret(newdim * hf.size());
    generic_compress(in, ret, newdim, hf, applicator);
    return ret;
}
template<typename C, typename C2, typename Hasher>
auto &cs_compress(const C &in, C2 &ret, size_t newdim, const Hasher &hf) {
    return generic_compress(in, ret, newdim, hf, CountApplicator());
}
template<typename C, typename Hasher=KWiseHasherSet<4>>
auto cs_compress(const C &in, size_t newdim, const Hasher &hf) {
    return generic_compress(in, newdim, hf, CountApplicator());
}

// Note: cs_compress could be a special case of wz_compress with Samplingist always returning 1
// or make it always a function of the remainder.

template<typename C, typename C2, typename Hasher=KWiseHasherSet<4>, typename SamplingDist=std::exponential_distribution<double>,
         typename RNG=blaze::RNG>
auto &wz_compress(const C &in, C2 &out, size_t newdim, const Hasher &hf, double p) {
    return generic_compress(in, out, newdim, hf, WZApplicator(p));
#if 0
    //using FT = std::decay_t<decltype(*std::begin(in))>;
    std::fill(out.begin(), out.end(), static_cast<std::decay_t<decltype(out[0])>>(0));
    if(newdim > in.size()) throw 1;
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(newdim);
    SamplingDist gen(p);
    for(unsigned j = 0; j < ns; ++j) {
        for(unsigned i = 0; i < in.size(); ++i) {
            const auto v = in[i];
            auto hv = hf(i, j);
            auto dm = div.divmod(hv);
            auto ind = dm.rem * ns + j;
            RNG rng(dm.quot >> 1);
            const double mult = gen(rng) * (dm.quot & 1 ? 1: -1);
            out.operator[](ind) += v * mult;
        }
        // TODO: decompress.
        // Sample using the same seed, just multiply by inverse
    }
    return out;
#endif
}

template<typename FT, bool SO, typename C2, typename Hasher=KWiseHasherSet<4>, typename SamplingDist=std::exponential_distribution<double>,
         typename RNG=blaze::RNG>
auto &wz_compress(const blaze::CompressedVector<FT> &in, C2 &out, size_t newdim, const Hasher &hf, double p) {
    return generic_compress(in, out, newdim, hf, WZApplicator(p));
#if 0
    //using FT = std::decay_t<decltype(*std::begin(in))>;
    std::fill(out.begin(), out.end(), static_cast<std::decay_t<decltype(out[0])>>(0));
    if(newdim > in.size()) throw 1;
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(newdim);
    SamplingDist gen(p);
    for(const auto &p: in) {
        for(unsigned j = 0; j < ns; ++j) {
            const auto v = p.value();
            auto hv = hf(p.index(), j);
            auto dm = div.divmod(hv);
            auto ind = dm.rem * ns + j;
            RNG rng(dm.quot >> 1);
            const double mult = gen(rng) * (dm.quot & 1 ? 1: -1);
            out.operator[](ind) += v * mult;
        }
        // TODO: decompress.
        // Sample using the same seed, just multiply by inverse
    }
    return out;
#endif
}

template<typename C, typename C2=C, typename Hasher=KWiseHasherSet<4>, typename RNG=blaze::RNG>
auto wz_compress(const C &in, size_t newdim, const Hasher &hf, double p) {
    C2 ret(newdim * hf.size());
    wz_compress(in, ret, newdim, hf, p);
    return ret;
}

template<typename C, typename OutC, typename Hasher=KWiseHasherSet<4>,
         typename=std::enable_if_t<!std::is_arithmetic<OutC>::value>,
         typename RNG=blaze::RNG>
auto &wz_decompress(const C &in, const Hasher &hf, OutC &ret, double p) {
    // TODO: some kind of importance sampling, weighting larger items more.
    PREC_REQ(in.size() % hf.size() == 0, "in dimension must be divisible by hf count");
    size_t olddim = in.size() / hf.size();
    size_t newdim = ret.size();
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(olddim);
    std::exponential_distribution<double> gen(p);
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < newdim; ++i) {
        sketch::common::detail::tmpbuffer<float, 9> mem(hf.size());
        auto tmp = mem.get();
        for(unsigned j = 0; j < ns; ++j) {
            auto hv = hf(i, j);
            auto dm = div.divmod(hv);
            RNG rng(dm.quot >> 1);
            tmp[j] = in[dm.rem * ns + j] / gen(rng) * (dm.quot & 1 ? 1: -1);
        }
        common::sort::insertion_sort(tmp, tmp + hf.size());
        ret[i] = (tmp[ns >> 1] + tmp[(ns - 1) >> 1]) * .5;
    }
    return ret;
}

template<typename C, typename OutC=C, typename Hasher=KWiseHasherSet<4>, typename RNG=blaze::RNG>
auto wz_decompress(const C &in, size_t outdim, const Hasher &hs, double p) {
    OutC ret(outdim);
    wz_decompress(in, hs, ret, p);
    return ret;
}

template<typename C, typename OutC, typename Hasher=KWiseHasherSet<4>,
         typename=std::enable_if_t<!std::is_arithmetic<Hasher>::value>>
auto &cs_decompress(const C &in, const Hasher &hf, OutC &ret) {
    PREC_REQ(in.size() % hf.size() == 0, "in dimension must be divisible by hf count");
    size_t olddim = in.size() / hf.size();
    size_t newdim = ret.size();
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(olddim);
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < newdim; ++i) {
        sketch::common::detail::tmpbuffer<float, 9> mem(hf.size());
        auto tmp = mem.get();
        for(unsigned j = 0; j < ns; ++j) {
            auto hv = hf(i, j);
            tmp[j] = in[div.mod(hv >> 1) * ns + j] * (hv & 1 ? 1: -1);
        }
        common::sort::insertion_sort(tmp, tmp + hf.size());
        ret[i] = (tmp[ns >> 1] + tmp[(ns - 1) >> 1]) * .5;
    }
    return ret;
}

template<typename C, typename OutC=C, typename Hasher=KWiseHasherSet<4>>
auto cs_decompress(const C &in, size_t newdim, const Hasher &hf) {
    OutC ret(newdim);
    cs_decompress<C, OutC, Hasher>(in, hf, ret);
    return ret;
}

struct AbsMax {
    template<typename T>
    bool operator()(T x, T y) const {return std::abs(x) > std::abs(y);}
};

template<typename C, typename Functor=std::greater<void>>
auto top_indices_from_compressed(const C &in, size_t newdim, size_t olddim, const KWiseHasherSet<4> &hf, unsigned k) {
    //if(newdim < in.size()) throw 1;
    const size_t ns = hf.size();
    schism::Schismatic<uint32_t> div(olddim);
    using FT = std::decay_t<decltype(*std::begin(in))>;
    std::priority_queue<std::pair<FT, unsigned>, std::vector<std::pair<FT, unsigned>>, Functor> pq;
    OMP_PRAGMA("omp parallel for")
    for(size_t i = 0; i < newdim; ++i) {
        sketch::common::detail::tmpbuffer<float, 8> mem(hf.size());
        auto tmp = mem.get();
        for(unsigned j = 0; j < ns; ++j) {
            auto hv = hf(i, j);
            tmp[j] = in.operator[](div.mod(hv >> 1) * ns + j) * (hv & 1 ? 1: -1);
        }
        common::sort::insertion_sort(tmp, tmp + hf.size());
        std::pair<FT, unsigned> pair = std::make_pair((tmp[ns >> 1] + tmp[(ns - 1) >> 1]) * .5, unsigned(i));
        OMP_PRAGMA("omp critical")
        {
            pq.push(pair);
            if(pq.size() > k) pq.pop();
        }
    }
    std::pair<std::vector<FT>, std::vector<unsigned>> ret;
    ret.first.reserve(k);
    ret.second.reserve(k);
    for(unsigned i = 0; i < k; ++i) {
        const auto &r(pq.top());
        ret.first.push_back(r.first);
        ret.second.push_back(r.second);
        pq.pop();
    }
    return ret;
}

template<typename HasherSetType=XORSeedHasherSet<>>
struct SketchApplicator {
protected:
    size_t in_, out_;
    double p_;
    const Transform tx_;
    HasherSetType hs_;
public:
    SketchApplicator(size_t indim, size_t outdim, uint64_t seed=13, double p=1., Transform tx=CountSketch):
        in_(indim), out_(outdim), p_(p), tx_(tx), hs_(seed)
    {
        PREC_REQ(tx == CountSketch || tx == WoodruffZhang, "Unsupported");
    }
    template<typename C, typename C2>
    auto &compress(const C &in, C2 &out) const {
        assert(in_ == in.size());
        switch(tx_) {
            case CountSketch:
                return cs_compress(in, out, out_, hs_);
            case WoodruffZhang:
                return wz_compress(in, out, out_, hs_, p_);
            default: __builtin_unreachable();
        }
    }
    template<typename C, typename C2>
    auto &decompress(const C &in, C2 &out) const {
        assert(out_ * hs_.size() == in.size());
        switch(tx_) {
            case CountSketch:
                return cs_decompress(in, hs_, out);
            case WoodruffZhang:
                return wz_decompress(in, hs_, out, p_);
            default: __builtin_unreachable();
        }
    }
    template<typename C>
    auto decompress(const C &in) const {
        assert(out_ * hs_.size() == in.size());
        switch(tx_) {
            case CountSketch:
                return cs_decompress(in, out_, hs_);
            case WoodruffZhang:
                return wz_decompress(in, out_, hs_, p_);
            default: __builtin_unreachable();
        }
    }
    template<typename C>
    auto compress(const C &in) const {
        assert(in_ == in.size());
        switch(tx_) {
            case CountSketch:
                return cs_compress(in, out_, hs_);
            case WoodruffZhang:
                return wz_compress(in, out_, hs_, p_);
            default: __builtin_unreachable();
        }
    }
};

} // rnla

} // sketch
