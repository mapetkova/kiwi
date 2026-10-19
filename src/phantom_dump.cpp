#include "kiwi/phantom_dump.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace kiwi {

namespace {

constexpr int32_t kMagic1 = 60769, kMagic2 = 60878, kMagic3 = 690706;
constexpr int kTagLen = 16;
constexpr int kNumKinds = 8;

// Reads a T from `p`, byte-swapping first if the file's byte order differs from the host's.
template <class T>
T load_value(const char* p, bool swap) {
    unsigned char buf[sizeof(T)];
    std::memcpy(buf, p, sizeof(T));
    if (swap) std::reverse(buf, buf + sizeof(T));
    T v;
    std::memcpy(&v, buf, sizeof(T));
    return v;
}

// One of the eight data types of the format (default int, int*1, int*2, int*4, int*8, default real,
// real*4, real*8), reduced to its byte size and whether it is floating point.
struct DataKind {
    int size = 0;
    bool is_float = false;
};

double decode(const char* p, const DataKind& k, bool swap) {
    if (k.is_float) return k.size == 4 ? static_cast<double>(load_value<float>(p, swap)) : load_value<double>(p, swap);
    switch (k.size) {
        case 1: return static_cast<double>(load_value<int8_t>(p, swap));
        case 2: return static_cast<double>(load_value<int16_t>(p, swap));
        case 4: return static_cast<double>(load_value<int32_t>(p, swap));
        default: return static_cast<double>(load_value<int64_t>(p, swap));
    }
}

// Fortran unformatted sequential access: every record is [4-byte length][payload][4-byte length].
class RecordReader {
public:
    explicit RecordReader(const std::string& path) : path_(path), f_(path, std::ios::binary) {}
    bool ok() const { return static_cast<bool>(f_); }
    bool swap() const { return swap_; }

    // Decodes the first record as the Phantom capture pattern (i1 r1 i2 iversion i3) and, if it matches,
    // fixes the file's byte order and the default int/real widths. Returns false if it is not a Phantom
    // dump. The record length is 4*int_bytes + real_bytes, which is also how the byte order is found:
    // only one of the two orders turns the first four bytes into one of the four possible lengths.
    bool read_capture(int& int_bytes, int& real_bytes, int& iversion) {
        char b[4];
        f_.read(b, 4);
        if (f_.gcount() != 4) return false;
        for (bool swap : {false, true}) {
            const int32_t n = load_value<int32_t>(b, swap);
            int ib = 0, rb = 0;
            if (n == 20) { ib = 4; rb = 4; }
            else if (n == 24) { ib = 4; rb = 8; }
            else if (n == 36) { ib = 8; rb = 4; }
            else if (n == 40) { ib = 8; rb = 8; }
            else continue;
            std::vector<char> body(static_cast<size_t>(n) + 4);
            f_.read(body.data(), static_cast<std::streamsize>(body.size()));
            if (static_cast<size_t>(f_.gcount()) != body.size()) return false;
            if (load_value<int32_t>(body.data() + n, swap) != n) return false;
            const char* p = body.data();
            auto get_int = [&](int off) { return ib == 4 ? static_cast<int64_t>(load_value<int32_t>(p + off, swap))
                                                          : load_value<int64_t>(p + off, swap); };
            const int64_t i1 = get_int(0);
            const double r1 = rb == 4 ? static_cast<double>(load_value<float>(p + ib, swap)) : load_value<double>(p + ib, swap);
            const int64_t i2 = get_int(ib + rb);
            const int64_t iv = get_int(2 * ib + rb);
            const int64_t i3 = get_int(3 * ib + rb);
            if (i1 != kMagic1 || i2 != kMagic2 || i3 != kMagic3 || r1 != static_cast<double>(kMagic2)) return false;
            swap_ = swap;
            int_bytes = ib;
            real_bytes = rb;
            iversion = static_cast<int>(iv);
            return true;
        }
        return false;
    }

    // Payload of the next record. Throws on truncation or mismatched markers.
    std::vector<char> record(const std::string& what) {
        const int32_t n = read_marker(what);
        std::vector<char> payload(static_cast<size_t>(n));
        if (n > 0) read_bytes(payload.data(), static_cast<size_t>(n), what);
        const int32_t m = read_marker(what);
        if (m != n) fail("record markers disagree (" + std::to_string(n) + " vs " + std::to_string(m) + ")", what);
        return payload;
    }

    // Skips the next record (arrays that are not needed are never decoded).
    void skip(const std::string& what) {
        const int32_t n = read_marker(what);
        f_.seekg(static_cast<std::streamoff>(n), std::ios::cur);
        const int32_t m = read_marker(what);
        if (m != n) fail("record markers disagree while skipping", what);
    }

    [[noreturn]] void fail(const std::string& why, const std::string& what) const {
        throw std::runtime_error("phantom-binary snapshot: " + why + " while reading " + what + ": " + path_);
    }

private:
    int32_t read_marker(const std::string& what) {
        char b[4];
        read_bytes(b, 4, what);
        const int32_t n = load_value<int32_t>(b, swap_);
        if (n < 0) fail("negative record length (records over 2 GB are not supported)", what);
        return n;
    }
    void read_bytes(char* dst, size_t n, const std::string& what) {
        f_.read(dst, static_cast<std::streamsize>(n));
        if (static_cast<size_t>(f_.gcount()) != n) fail("unexpected end of file", what);
    }

    std::string path_;
    std::ifstream f_;
    bool swap_ = false;
};

std::string trim_tag(const char* p, size_t n) {
    std::string s(p, n);
    const size_t a = s.find_first_not_of(' ');
    if (a == std::string::npos) return std::string();
    const size_t b = s.find_last_not_of(' ');
    return s.substr(a, b - a + 1);
}

} // namespace

bool is_phantom_dump(const std::string& path) {
    try {
        RecordReader r(path);
        if (!r.ok()) return false;
        int ib, rb, iv;
        return r.read_capture(ib, rb, iv);
    } catch (...) {
        return false;
    }
}

SPHSnapshot load_phantom_dump(const std::string& path, const PhantomUnits& units) {
    RecordReader r(path);
    if (!r.ok()) throw std::runtime_error("could not open SPH snapshot file: " + path);

    int int_bytes = 0, real_bytes = 0, iversion = 0;
    if (!r.read_capture(int_bytes, real_bytes, iversion))
        throw std::runtime_error("phantom-binary snapshot: the file does not start with the Phantom capture "
                                 "pattern (is it a Phantom dump?): " + path);
    const bool swap = r.swap();

    // 100-character identifier: first letter F (full dump) or S (small dump), second T (tagged format).
    const std::vector<char> ident = r.record("the file identifier");
    const std::string ident_str(ident.begin(), ident.end());
    if (ident.size() < 2 || (ident[0] != 'F' && ident[0] != 'S') || ident[1] != 'T')
        r.fail("only tagged full dumps (file identifier 'FT:...') are supported, got '" + ident_str.substr(0, 20) + "'",
               "the file identifier");

    const DataKind kinds[kNumKinds] = {{int_bytes, false}, {1, false}, {2, false}, {4, false},
                                       {8, false},         {real_bytes, true}, {4, true}, {8, true}};

    // ---- global header: eight typed sections of (tag, value) pairs; a tag may repeat (one entry per type)
    std::map<std::string, std::vector<double>> header;
    for (int k = 0; k < kNumKinds; ++k) {
        const std::vector<char> cnt = r.record("a header section size");
        if (cnt.size() != 4) r.fail("malformed header section size", "the global header");
        const int32_t nvars = load_value<int32_t>(cnt.data(), swap);
        if (nvars <= 0) continue;
        const std::vector<char> tags = r.record("header tags");
        const std::vector<char> vals = r.record("header values");
        if (tags.size() != static_cast<size_t>(nvars) * kTagLen ||
            vals.size() != static_cast<size_t>(nvars) * static_cast<size_t>(kinds[k].size))
            r.fail("header section " + std::to_string(k) + " has inconsistent sizes", "the global header");
        for (int i = 0; i < nvars; ++i)
            header[trim_tag(tags.data() + i * kTagLen, kTagLen)].push_back(
                decode(vals.data() + static_cast<size_t>(i) * kinds[k].size, kinds[k], swap));
    }
    auto need = [&](const char* key) -> const std::vector<double>& {
        auto it = header.find(key);
        if (it == header.end() || it->second.empty())
            r.fail(std::string("required header entry '") + key + "' is missing", "the global header");
        return it->second;
    };
    const std::vector<double>& massoftype = need("massoftype");
    const double hfact = need("hfact")[0];
    const double udist = need("udist")[0], umass = need("umass")[0], utime = need("utime")[0];
    if (hfact <= 0.0 || udist <= 0.0 || umass <= 0.0 || utime <= 0.0)
        r.fail("header hfact/udist/umass/utime must be > 0", "the global header");
    const int hdr_blocks = header.count("nblocks") ? std::max(1, static_cast<int>(header["nblocks"][0])) : 1;

    // ---- array block headers
    const std::vector<char> nb = r.record("the number of array blocks");
    if (nb.size() != 4) r.fail("malformed block count", "the array section");
    const int nblocks = load_value<int32_t>(nb.data(), swap);
    if (nblocks <= 0 || nblocks % hdr_blocks != 0)
        r.fail("array block count " + std::to_string(nblocks) + " is not a multiple of the header's nblocks " +
                   std::to_string(hdr_blocks), "the array section");
    struct BlockHeader { int64_t n; int32_t nums[kNumKinds]; };
    std::vector<BlockHeader> blocks(static_cast<size_t>(nblocks));
    for (auto& b : blocks) {
        const std::vector<char> h = r.record("a block header");
        if (h.size() != 8 + 4 * kNumKinds) r.fail("malformed block header", "the array section");
        b.n = load_value<int64_t>(h.data(), swap);
        for (int k = 0; k < kNumKinds; ++k) b.nums[k] = load_value<int32_t>(h.data() + 8 + 4 * k, swap);
        if (b.n < 0) r.fail("negative particle count", "the array section");
    }

    // Particle types: with several types present, the itype array says which is which.
    const int ntypes = header.count("ntypes") ? static_cast<int>(header["ntypes"][0]) : static_cast<int>(massoftype.size());
    bool several_types = false;
    if (header.count("npartoftype")) {
        const auto& npt = header["npartoftype"];
        for (int t = 1; t < ntypes && t < static_cast<int>(npt.size()); ++t)
            if (npt[t] > 0.0) several_types = true;
    }

    const std::set<std::string> wanted = {"x", "y", "z", "vx", "vy", "vz", "h", "itype"};
    const double lf = udist / units.length_cm;                 // code length -> output length
    const double mf = umass / units.mass_g;                    // code mass   -> output mass
    const double vf = (udist / utime) / units.velocity_cms;    // code velocity -> output velocity
    const double rho_cgs = umass / (udist * udist * udist);    // code density -> g/cm^3
    const double m_gas_code = massoftype[0];

    SPHSnapshot snap;
    long long n_seen = 0, n_dead = 0, n_other_type = 0;
    const int nsub = nblocks / hdr_blocks;   // blocks per MPI task: [particles, sinks] normally
    for (int bi = 0; bi < nblocks; ++bi) {
        const BlockHeader& bh = blocks[static_cast<size_t>(bi)];
        const bool is_sink_block = (bi % nsub) == 1;
        std::map<std::string, std::vector<double>> arr;
        for (int k = 0; k < kNumKinds; ++k) {
            for (int a = 0; a < bh.nums[k]; ++a) {
                const std::vector<char> tag = r.record("an array tag");
                if (tag.size() != static_cast<size_t>(kTagLen)) r.fail("malformed array tag", "the array section");
                const std::string name = trim_tag(tag.data(), kTagLen);
                if (is_sink_block || !wanted.count(name) || arr.count(name)) {
                    r.skip("array '" + name + "'");
                    continue;
                }
                const std::vector<char> data = r.record("array '" + name + "'");
                if (data.size() != static_cast<size_t>(bh.n) * static_cast<size_t>(kinds[k].size))
                    r.fail("array '" + name + "' has " + std::to_string(data.size()) + " bytes, expected " +
                               std::to_string(bh.n * kinds[k].size), "the array section");
                std::vector<double>& v = arr[name];
                v.resize(static_cast<size_t>(bh.n));
                for (int64_t i = 0; i < bh.n; ++i)
                    v[static_cast<size_t>(i)] = decode(data.data() + i * kinds[k].size, kinds[k], swap);
            }
        }
        if (is_sink_block) continue;
        if (bh.n == 0) continue;
        for (const char* name : {"x", "y", "z", "vx", "vy", "vz", "h"})
            if (!arr.count(name))
                r.fail(std::string("array '") + name + "' not found in a particle block (small dumps do not "
                       "carry every array KIWI needs)", "the array section");
        const bool have_itype = arr.count("itype") != 0;
        if (!have_itype && several_types)
            r.fail("the header lists particles of more than one type but the dump has no itype array",
                   "the array section");
        for (int64_t i = 0; i < bh.n; ++i) {
            ++n_seen;
            const size_t u = static_cast<size_t>(i);
            const double h = arr["h"][u];
            if (!(h > 0.0)) { ++n_dead; continue; }   // h <= 0: accreted onto a sink / dead
            if (have_itype && std::lround(arr["itype"][u]) != 1) { ++n_other_type; continue; }
            snap.id.push_back(static_cast<int>(snap.id.size()) + 1);
            snap.x.push_back(arr["x"][u] * lf);
            snap.y.push_back(arr["y"][u] * lf);
            snap.z.push_back(arr["z"][u] * lf);
            snap.vx.push_back(arr["vx"][u] * vf);
            snap.vy.push_back(arr["vy"][u] * vf);
            snap.vz.push_back(arr["vz"][u] * vf);
            snap.m.push_back(m_gas_code * mf);
            snap.h.push_back(h * lf);
            const double q = hfact / h;
            snap.rho.push_back(m_gas_code * q * q * q * rho_cgs);   // g/cm^3, like the splash_ascii density column
        }
    }
    if (snap.size() == 0) r.fail("no gas particles found", "the array section");
    std::cerr << "[kiwi] phantom-binary snapshot: kept " << snap.size() << " of " << n_seen
              << " particle(s) (dropped " << n_dead << " dead with h <= 0 and " << n_other_type
              << " of another type; sink particles are ignored); gas mass " << m_gas_code * mf
              << " each, density from the header's hfact = " << hfact << ": " << path << "\n";
    return snap;
}

} // namespace kiwi
