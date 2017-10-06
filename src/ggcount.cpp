/*
 * Computing the Number of Paths in a Grid Graph
 * Hiroaki Iwashita <iwashita@erato.ist.hokudai.ac.jp>
 * Copyright (c) 2013 ERATO MINATO Project
 * $Id$
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef BIGNUM_SIZE
#define BIGNUM_SIZE 14
#endif

typedef uint64_t Code;
typedef uint64_t MateID;

constexpr int MATE_SIZE = 4 * sizeof(MateID);

enum {
    NONE, DEFAULT, VERBOSE
} msg;

int rows = 0;
int cols = 0;
uint64_t modulus = 0;
bool modnumMode = false;
bool cycleMode = false;
bool hamiltMode = false;
double startTime = 0;

void usage(char const* cmd) {
    std::cerr << "Usage: " << cmd
            << " [<option>...] <columns> [<rows>] [%<modulus>]\n";
    std::cerr << "Options\n";
    std::cerr << "  -c : Count cycles instead of paths\n";
    std::cerr << "  -h : Count Hamiltonian paths/cycles\n";
    std::cerr << "  -v : Print verbose messages\n";
    std::cerr << "  -q : Work quietly\n";
}

int bitcount(uint32_t n) {
    n = ((n & 0xaaaaaaaa) >> 1) + (n & 0x55555555);
    n = ((n & 0xcccccccc) >> 2) + (n & 0x33333333);
    n = ((n & 0xf0f0f0f0) >> 4) + (n & 0x0f0f0f0f);
    n = ((n & 0xff00ff00) >> 8) + (n & 0x00ff00ff);
    return ((n & 0xffff0000) >> 16) + (n & 0x0000ffff);
}

template<int size>
class Bignum {
    uint32_t val[size];

public:
    Bignum() {
    }

    Bignum(uint32_t i) {
        val[0] = i;
        for (int i = 1; i < size; ++i) {
            val[i] = 0;
        }
    }

    Bignum& operator=(uint32_t i) {
        val[0] = i;
        for (int i = 1; i < size; ++i) {
            val[i] = 0;
        }
        return *this;
    }

    bool operator==(Bignum const& o) const {
        for (int i = 0; i < size; ++i) {
            if (val[i] != o.val[i]) return false;
        }
        return true;
    }

    bool operator!=(Bignum const& o) const {
        return !operator==(o);
    }

    bool operator==(uint32_t i) const {
        if (val[0] != i) return false;
        for (int i = 1; i < size; ++i) {
            if (val[i] != 0) return false;
        }
        return true;
    }

    bool operator!=(uint32_t i) const {
        return !operator==(i);
    }

    void operator+=(Bignum const& o) {
        uint64_t x = 0;
        for (int i = 0; i < size; ++i) {
            x += val[i];
            x += o.val[i];
            val[i] = x;
            x >>= 32;
        }
        if (x != 0) throw std::runtime_error(
                "Bignum<" + std::to_string(size) + "> overflow!");
    }

    uint32_t divide(uint32_t n) {
        uint64_t r = 0;
        for (int i = size - 1; i >= 0; --i) {
            r = (r << 32) + val[i];
            lldiv_t d = lldiv(r, n);
            val[i] = d.quot;
            r = d.rem;
        }
        return r;
    }

    friend std::ostream& operator<<(std::ostream& os, Bignum const& o) {
        Bignum n = o;
        n.printHelper(os);
        return os;
    }

private:
    void printHelper(std::ostream& os) {
        uint32_t r = divide(10);
        if (*this != 0) printHelper(os);
        os << r;
    }
};

template<typename T>
class Modnum {
    T val;

public:
    Modnum() {
    }

    Modnum(T const& i)
            : val(i) {
        assert(i < modulus);
    }

    Modnum& operator=(T const& i) {
        assert(i < modulus);
        val = i;
        return *this;
    }

    void operator+=(Modnum const& o) {
        assert(0 <= val && val < modulus);
        assert(0 <= o.val && o.val < modulus);
        val += o.val;
        if (val >= modulus || val < o.val) val -= modulus;
        assert(0 <= val && val < modulus);
    }

    operator T() const {
        return val;
    }

    friend std::ostream& operator<<(std::ostream& os, Modnum const& o) {
        return os << o.val + 0;
    }
};

enum MateValue {
    N = 0, R = 1, L = 2, X = 3
};

enum MateValuePair {
    NN = 0x0, NR = 0x1, NL = 0x2, NX = 0x3, //
    RN = 0x4, RR = 0x5, RL = 0x6, RX = 0x7, //
    LN = 0x8, LR = 0x9, LL = 0xa, LX = 0xb, //
    XN = 0xc, XR = 0xd, XL = 0xe, XX = 0xf, //
};

std::ostream& operator<<(std::ostream& os, MateValue v) {
    static char const* tbl = ".)(X";
    return os << tbl[v];
}

std::ostream& operator<<(std::ostream& os, MateValuePair w) {
    static char const* tbl = ".)(X";
    return os << tbl[w >> 2] << tbl[w & 3];
}

class Mate {
    MateID value;

    static MateID mask(int k, MateValue u) {
        return MateID(u) << (2 * k);
    }

    static MateID mask(int k, MateValuePair w) {
        return MateID(w) << (2 * (k - 1));
    }

public:
    Mate()
            : value(0) {
    }

    Mate(MateID id)
            : value(id) {
    }

    Mate(int k, MateValue u)
            : value(mask(k, u)) {
    }

    Mate(int k, MateValuePair w)
            : value(mask(k, w)) {
    }

    MateID id() const {
        return value;
    }

    MateValue get(int k) const {
        return MateValue((value >> (2 * k)) & 0x3);
    }

    void set(int k, MateValue v) {
        value = (value & ~mask(k, X)) | mask(k, v);
    }

    MateValuePair getPair(int k) const {
        return MateValuePair((value >> (2 * (k - 1))) & 0xf);
    }

    void setPair(int k, MateValuePair w) {
        value = (value & ~mask(k, XX)) | mask(k, w);
    }

    Mate shrink(int k) const {
        MateID m = (MateID(1) << (2 * k)) - 1;
        MateID h = value & ~m;
        MateID l = value & m;
        return (h >> 2) | l;
    }

    Mate getRight(int k) const {
        return value & ((MateID(1) << (k * 2)) - 1);
    }

    Mate shiftLeft(int k) const {
        return value << (k * 2);
    }

    Mate shiftRight(int k) const {
        return value >> (k * 2);
    }

    Mate operator|(Mate const& o) const {
        return value | o.value;
    }

    Mate const& operator++() {
        ++value;
        if ((value & 3) == 3) {
            MateID v1 = 1;
            MateID v3 = 3;
            do {
                value += v1;
                v1 <<= 2;
                v3 <<= 2;
            } while (v3 != 0 && (value & v3) == v3);
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, Mate const& mate) {
        os << "[";
        for (int i = MATE_SIZE - 1; i >= 0; --i) {
            os << mate.get(i);
        }
        return os << "]";
    }
};

struct CodeTable {
    Code base;
    Code size;
    Mate mateL;
    Mate const* mateR;
};

class MateCodec {
    struct Block {
        Code size;
        Mate* mate;

        ~Block() {
            delete[] mate;
        }
    };

    int wL;
    int wR;
    int hL;
    int hR;
    int hC;

    MateID mateSizeL_;
    MateID mateSizeR_;
    Code codeSize_;
    Code codeSizeL_;
    Code* codeLR_;
    Code* codeL_;
    Code* codeR_;
    Block* rightBlock_;
    CodeTable* codeTable_;

public:
    MateCodec(int width, int rightWidth, int leftHeight, int rightHeight) {

        assert(width >= 1);
        wR = (rightWidth < 1) ? 1 : (width < rightWidth) ? width : rightWidth;
        wL = width - wR;
        hL = leftHeight;
        hR = rightHeight;
        hC = std::min(hL + wL, hR + wR);

        mateSizeL_ = (wL >= 1) ? Mate(wL - 1, X).id() : 1;
        mateSizeR_ = Mate(wR - 1, X).id();
        codeSize_ = 0;
        codeSizeL_ = 0;
        codeLR_ = new Code[mateSizeL_]();
        codeL_ = new Code[mateSizeL_]();
        codeR_ = new Code[mateSizeR_]();

        if (msg == VERBOSE) std::cerr << "Allocated (2 × " << mateSizeL_
                << " + " << mateSizeR_ << ") × " << sizeof(Code) << " = "
                << (mateSizeL_ * 2 + mateSizeR_) * sizeof(Code)
                << " bytes for code tables.\n";

        rightBlock_ = new Block[hC + 1]();

        for (int h = 0; h <= hC; ++h) {
            Code n = motzkin(wR, h, hR);
            rightBlock_[h].mate = new Mate[n]();

            if (msg == VERBOSE) std::cerr << "Allocated " << n << " × "
                    << sizeof(Mate) << " = " << n * sizeof(Mate)
                    << " bytes for right state list #" << h << ".\n";

            fillCodeR(h, wR, h, hR, Mate());

            assert(rightBlock_[h].size == n);
        }

        {
            Code n = 0;
            for (int h = 0; h <= hC; ++h) {
                n += motzkin(wL, hL, h);
            }
            codeTable_ = new CodeTable[n]();

            if (msg == VERBOSE) std::cerr << "Allocated " << n << " × "
                    << sizeof(CodeTable) << " = " << n * sizeof(CodeTable)
                    << " bytes for left state list.\n";

            fillCodeL(wL, hL, Mate());

            assert(codeSizeL_ == n);
        }

        //std::cerr << "\n" << *this;
    }

    ~MateCodec() {
        delete[] codeLR_;
        delete[] codeL_;
        delete[] codeR_;
        delete[] rightBlock_;
        delete[] codeTable_;
    }

    int leftWidth() const {
        return wL;
    }

    int rightWidth() const {
        return wR;
    }

    int leftHeight() const {
        return hL;
    }

    int rightHeight() const {
        return hR;
    }

    int centerHeight() const {
        return hC;
    }

    Code codeSize() const {
        return codeSize_;
    }

    Code codeSizeL() const {
        return codeSizeL_;
    }

    CodeTable const& codeTable(int i) const {
        return codeTable_[i];
    }

    Code encode(Mate const& mate) const {
        MateID mL = mate.shiftRight(wR).id();
        MateID mR = mate.getRight(wR).id();
        assert(mL < mateSizeL_);
        assert(mR < mateSizeR_);
        assert(mL == 0 || codeL_[mL] > 0);
        return codeLR_[mL] + codeR_[mR];
    }

    Code encodeL(Mate const& mateL) const {
        MateID mL = mateL.id();
        assert(mL < mateSizeL_);
        assert(mL == 0 || codeL_[mL] > 0);
        return codeL_[mL];
    }

    friend std::ostream& operator<<(std::ostream& os, MateCodec const& o) {
        for (Code i = 0; i < o.codeSizeL_; ++i) {
            Code c = o.codeTable_[i].base;
            Code n = o.codeTable_[i].size;
            Mate mL = o.codeTable_[i].mateL;
            Mate const* mR = o.codeTable_[i].mateR;
            if (mR == 0) continue;

            os << "block " << i << ": " << c << ".." << c + n - 1 << "\n";
            for (Code j = 0; j < n; ++j) {
                os << "  " << (mL | *mR++) << "\n";
            }
        }
        return os;
    }

private:
    static size_t motzkin(int w, int h1, int h2) {
        static std::vector<std::vector<std::vector<size_t>>> cache;

        if (w == 0 && h1 == h2) return 1;
        if (w <= 0 || h1 < 0 || h2 < 0) return 0;
        if (w < std::abs(h1 - h2)) return 0;

        if (cache.size() <= size_t(w)) cache.resize(w + 1);
        if (cache[w].size() <= size_t(h1)) cache[w].resize(h1 + 1);
        if (cache[w][h1].size() <= size_t(h2)) cache[w][h1].resize(h2 + 1);

        size_t& m = cache[w][h1][h2];
        if (m == 0) {
            m = motzkin(w - 1, h1 + 1, h2) + motzkin(w - 1, h1, h2)
                    + motzkin(w - 1, h1 - 1, h2);
        }
        return m;
    }

    void fillCodeR(int h, int w, int h1, int h2, Mate mR) {
        if (w == 0 && h1 == h2) {
            Block& rt = rightBlock_[h];
            codeR_[mR.id()] = rt.size;
            rt.mate[rt.size++] = mR;
            return;
        }
        if (w <= 0 || h1 < 0 || h2 < 0) return;
        if (w < std::abs(h1 - h2)) return;

        --w;
        fillCodeR(h, w, h1, h2, mR);
        mR.set(w, R);
        fillCodeR(h, w, h1 - 1, h2, mR);
        mR.set(w, L);
        fillCodeR(h, w, h1 + 1, h2, mR);
    }

    void fillCodeL(int w, int h, Mate mL) {
        if (w == 0 && 0 <= h && h <= hC) {
            Block& rt = rightBlock_[h];
            CodeTable& ct = codeTable_[codeSizeL_];
            ct.base = codeLR_[mL.id()] = codeSize_;
            ct.size = rt.size;
            ct.mateL = mL.shiftLeft(wR);
            ct.mateR = rt.mate;
            codeSize_ += rt.size;
            codeL_[mL.id()] = codeSizeL_++;
            return;
        }
        if (w <= 0 || h < 0 || w < h - hC) return;

        --w;
        fillCodeL(w, h, mL);
        mL.set(w, R);
        fillCodeL(w, h - 1, mL);
        mL.set(w, L);
        fillCodeL(w, h + 1, mL);
    }
};

template<typename Number>
class PathCounter {
    int const rows;
    int const cols;
    bool const cycleMode;
    bool const hamiltMode;
    MateCodec mc;
    MateCodec wc;

    Number* value;
    Number* deferred;

    Number total;

    int groupWidth;
    int numGroups;
    int* groups;

public:
    PathCounter(int rows, int cols, bool cycle, bool hamilt)
            : rows(rows), cols(cols), cycleMode(cycle), hamiltMode(hamilt),
              mc(cols, (cols + 1) / 2, cycle ? 0 : 1, 0),
              wc(cols - 1, cols / 2, cycle ? 0 : 1, 0) {

        value = new Number[mc.codeSize()];
        deferred = new Number[wc.codeSize()];

        if (mc.leftWidth() >= 3) {
            groupWidth = mc.leftWidth() - 2;
            numGroups = 1 << groupWidth;
            groups = new int[numGroups];

            if (msg == VERBOSE) std::cerr << "Allocated " << numGroups << " × "
                    << sizeof(int) << " = " << numGroups * sizeof(int)
                    << " bytes for group list.\n";
        }
        else {
            groupWidth = 0;
            numGroups = 0;
            groups = 0;
        }

        for (int i = 0; i < numGroups; ++i) {
            groups[i] = i;
        }

        std::sort(groups, groups + numGroups, [](int a, int b) {
            return bitcount(a) > bitcount(b);
        });

        if (msg) std::cerr << "Allocated (" << mc.codeSize() << " + "
                << wc.codeSize() << ") × " << sizeof(Number) << " = "
                << (mc.codeSize() + wc.codeSize()) * sizeof(Number)
                << " bytes for numbers.\n";
    }

    ~PathCounter() {
        delete[] value;
        delete[] deferred;
        delete[] groups;
    }

private:
    void update(int j, bool final) {
        int p = cols - j - 1; // bit position
        assert(1 <= p && p < cols);

        if (cycleMode && (!hamiltMode || final)) {
            total += value[mc.encode(Mate(p, LR))];
        }

        if (groups) {
            int ungroupPos =
                    (p - mc.rightWidth() > 1) ? p - mc.rightWidth() : 1;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (int i = 0; i < numGroups; ++i) {
                updateGroup(p, ungroupPos, groups[i], 0, mc.leftWidth(),
                        mc.leftHeight());
            }
        }
        else {
            for (Code i = 0; i < mc.codeSizeL(); ++i) {
                updateBlock(p, mc.codeTable(i));
            }
        }
    }

    void updateGroup(int p, int ungroupPos, int group, Mate mL, int w,
            int h) const {
        if (w == 0 && 0 <= h && h <= mc.centerHeight()) {
            Code iL = mc.encodeL(mL);
            updateBlock(p, mc.codeTable(iL));
            return;
        }
        if (w <= 0 || h < 0 || w < h - mc.centerHeight()) return; // out of bounds

        --w;
        if (w == ungroupPos || w == ungroupPos - 1) {
            mL.set(w, N), updateGroup(p, ungroupPos, group, mL, w, h);
            mL.set(w, R), updateGroup(p, ungroupPos, group, mL, w, h - 1);
            mL.set(w, L), updateGroup(p, ungroupPos, group, mL, w, h + 1);
        }
        else {
            int k = (w > ungroupPos) ? w - 2 : w;
            if ((group >> k) & 1) {
                mL.set(w, R), updateGroup(p, ungroupPos, group, mL, w, h - 1);
                mL.set(w, L), updateGroup(p, ungroupPos, group, mL, w, h + 1);
            }
            else {
                mL.set(w, N), updateGroup(p, ungroupPos, group, mL, w, h);
            }
        }
    }

    void updateBlock(int p, CodeTable const& block) const {
        for (Code i = 0; i < block.size; ++i) {
            Mate mate = block.mateL | block.mateR[i];
            Number& c = value[block.base + i];
            MateValuePair w = mate.getPair(p);

            switch (w) {
            case NN: // [..] → [()]
                {
                    Number& d = deferred[wc.encode(mate.shrink(p))];

                    if (c != 0) {
                        mate.setPair(p, LR);
                        value[mc.encode(mate)] += c;
                    }

                    if (hamiltMode) { // Can't leave untouched
                        c = d;
                    }
                    else {
                        c += d;
                    }

                    d = 0;
                }
                break;
            case NL: // [.(] ↔ [(.]
            case NR: // [.)] ↔ [).]
                {
                    Number& d = deferred[wc.encode(mate.shrink(p))];

                    mate.setPair(p, (w == NL) ? LN : RN);

                    Number& cc = value[mc.encode(mate)];

                    if (hamiltMode) { // Can't leave c untouched
                        if (p == 2 && mate.get(0) == N) {
                            // Can't take c
                            c = cc;
                            c += d;
                            d = 0;
                        }
                        else if (p == 1) {
                            // Can't leave cc untouched
                            d += cc;
                            cc = c;
                            c = d;
                            d = 0;
                        }
                        else {
                            Number tmp = c;
                            c = cc;
                            c += d;
                            d = tmp;
                        }
                    }
                    else {
                        if (p == 1) {
                            d += cc;
                            cc += c;
                            c += d;
                            d = 0;
                        }
                        else {
                            Number tmp = c;
                            c += cc;
                            c += d;
                            d = tmp;
                        }
                    }
                }
                break;
            case LL: // [((---)] → [..---(]
                {
                    mate.setPair(p, NN);

                    int q = p - 1;
                    int s = 1;
                    while (s > 0) {
                        --q;
                        assert(q >= 0);
                        switch (mate.get(q)) {
                        case L:
                            ++s;
                            break;
                        case R:
                            --s;
                            break;
                        default:
                            break;
                        }
                    }
                    mate.set(q, L);

                    if (p == 1) {
                        value[mc.encode(mate)] += c;
                    }
                    else {
                        deferred[wc.encode(mate.shrink(p - 1))] += c;
                    }
                }
                break;
            case RR: // [(---))] → [)---..]
                if (!(hamiltMode && p == 2 && mate.get(0) == N)) {
                    mate.setPair(p, NN);

                    int q = p;
                    int s = 1;
                    while (s > 0) {
                        ++q;
                        assert(q < cols);
                        switch (mate.get(q)) {
                        case L:
                            --s;
                            break;
                        case R:
                            ++s;
                            break;
                        default:
                            break;
                        }
                    }
                    mate.set(q, R);

                    if (p == 1) {
                        value[mc.encode(mate)] += c;
                    }
                    else {
                        deferred[wc.encode(mate.shrink(p - 1))] += c;
                    }
                }
                break;
            case RL: // [)(] → [..]
                if (!(hamiltMode && p == 2 && mate.get(0) == N)) {
                    mate.setPair(p, NN);
                    if (p == 1) {
                        value[mc.encode(mate)] += c;
                    }
                    else {
                        deferred[wc.encode(mate.shrink(p - 1))] += c;
                    }
                }
                break;
            default:
                break;
            }
        }
    }

public:
    Number count() {
        if (msg) std::cerr << "Counting" << (hamiltMode ? " Hamiltonian" : "")
                << (cycleMode ? " cycles" : " paths") << " in a " << cols << "x"
                << rows << " grid graph"
#ifdef _OPENMP
                << "; #thread = " << omp_get_max_threads()
#endif
                << ".\n";

        total = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (Code code = 0; code < mc.codeSize(); ++code) {
            value[code] = 0;
        }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (Code code = 0; code < wc.codeSize(); ++code) {
            deferred[code] = 0;
        }

        value[mc.encode(Mate(cols - 1, cycleMode ? N : R))] = 1;

        for (int i = 0; i < rows; ++i) {
            if (msg == DEFAULT) std::cerr
                    << (i == 0 ? "  ┌" : i < rows - 1 ? "  ├" : "  └");
            for (int j = 0; j < cols - 2; ++j) {
                update(j, false);
                if (msg == DEFAULT) {
                    std::cerr << "─"
                            << (i == 0 ? "┬" : i < rows - 1 ? "┼" : "┴");
                }
            }
            update(cols - 2, i == rows - 1);
            if (msg == DEFAULT) {
                std::cerr << "─" << (i == 0 ? "┐" : i < rows - 1 ? "┤" : "┘")
                        << "\n";
            }
        }

        if (!cycleMode) {
            total = value[mc.encode(Mate(0, R))];
        }

        return total;
    }
};

int main(int argc, char *argv[]) {
    msg = DEFAULT;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string s = argv[i];
            if (s[0] == '-') {
                for (size_t j = 1; j < s.length(); ++j) {
                    switch (s[j]) {
                    case 'c':
                        cycleMode = true;
                        break;
                    case 'h':
                        hamiltMode = true;
                        break;
                    case 'v':
                        msg = VERBOSE;
                        break;
                    case 'q':
                        msg = NONE;
                        break;
                    default:
                        throw std::exception();
                    }
                }
            }
            else if (s[0] == '%') {
                try {
                    modulus = std::stoull(s.substr(1));
                }
                catch (std::out_of_range e) {
                    modulus = 0;
                }
                modnumMode = true;
            }
            else if (cols == 0) {
                cols = std::stoi(s) + 1;
            }
            else if (rows == 0) {
                rows = std::stoi(s) + 1;
            }
            else {
                throw std::exception();
            }
        }

        if (rows == 0) rows = cols;
        if (rows <= 0 || cols <= 0) throw std::exception();
    }
    catch (std::exception& e) {
        usage(argv[0]);
        return 1;
    }

    if (rows < 1 || cols < 1) {
        std::cout << "0\n";
    }
    else if (rows == 1 || cols == 1) {
        std::cout << "1\n";
    }
    else if (!modnumMode) {
        PathCounter<Bignum<BIGNUM_SIZE>> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << "\n";
    }
    else if (modulus == 0) {
        std::string s = "0x1" + std::string(sizeof(uint64_t) * 2, '0');
        if (msg) std::cerr << ">>>>> MODULO " << s << " <<<<<\n";
        PathCounter<uint64_t> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod " << s << ")\n";
    }
    else if (modulus < 0x100) {
        if (msg) std::cerr << ">>>>> MODULO " << modulus << " <<<<<\n";
        PathCounter<Modnum<uint8_t>> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod " << modulus << ")\n";
    }
    else if (modulus == 0x100) {
        if (msg) std::cerr << ">>>>> MODULO 0x100 <<<<<\n";
        PathCounter<uint8_t> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() + 0 << " (mod 0x100)\n";
    }
    else if (modulus < 0x10000) {
        if (msg) std::cerr << ">>>>> MODULO " << modulus << " <<<<<\n";
        PathCounter<Modnum<uint16_t>> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod " << modulus << ")\n";
    }
    else if (modulus == 0x10000) {
        if (msg) std::cerr << ">>>>> MODULO 0x10000 <<<<<\n";
        PathCounter<uint16_t> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod 0x10000)\n";
    }
    else if (modulus < 0x100000000) {
        if (msg) std::cerr << ">>>>> MODULO " << modulus << " <<<<<\n";
        PathCounter<Modnum<uint32_t>> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod " << modulus << ")\n";
    }
    else if (modulus == 0x100000000) {
        if (msg) std::cerr << ">>>>> MODULO 0x100000000 <<<<<\n";
        PathCounter<uint32_t> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod 0x100000000)\n";
    }
    else {
        if (msg) std::cerr << ">>>>> MODULO " << modulus << " <<<<<\n";
        PathCounter<Modnum<uint64_t>> pc(rows, cols, cycleMode, hamiltMode);
        std::cout << pc.count() << " (mod " << modulus << ")\n";
    }

    return 0;
}
