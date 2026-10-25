#include "kiwi/snapshot.hpp"
#include "kiwi/phantom_dump.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace kiwi {

namespace {

bool is_comment_line(const std::string& line) {
    for (char c : line) {
        if (c == ' ' || c == '\t' || c == '\r') continue;
        return c == '#';
    }
    return false; // blank (or all-whitespace) line -- not a comment
}

bool is_blank_line(const std::string& line) {
    for (char c : line) {
        if (c != ' ' && c != '\t' && c != '\r') return false;
    }
    return true;
}

SPHSnapshot load_plain_ascii(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open SPH snapshot file: " + path);

    SPHSnapshot snap;
    // One particle per line, exactly 9 finite numbers. Parsed line by line (not as one token stream)
    // so that a malformed line is an error with its line number, rather than silently ending the
    // read there (dropping every later particle) or shifting the columns of the following rows.
    std::string line;
    long lineno = 0;
    int row = 0;
    while (std::getline(f, line)) {
        ++lineno;
        if (is_blank_line(line)) continue;
        double v[9];
        const char* p = line.c_str();
        for (int k = 0; k < 9; ++k) {
            char* end = nullptr;
            v[k] = std::strtod(p, &end);
            if (end == p)
                throw std::runtime_error(
                    "SPH snapshot file (plain_ascii format), line " + std::to_string(lineno) +
                    ": expected 9 numeric columns (x y z vx vy vz m h rho), could not read column " +
                    std::to_string(k + 1) + ": " + path);
            if (!std::isfinite(v[k]))
                throw std::runtime_error(
                    "SPH snapshot file (plain_ascii format), line " + std::to_string(lineno) +
                    ": column " + std::to_string(k + 1) + " is not finite: " + path);
            p = end;
        }
        while (*p == ' ' || *p == '\t' || *p == '\r') ++p;
        if (*p != '\0')
            throw std::runtime_error(
                "SPH snapshot file (plain_ascii format), line " + std::to_string(lineno) +
                ": more than 9 columns (expected x y z vx vy vz m h rho): " + path);
        ++row;
        snap.id.push_back(row);
        snap.x.push_back(v[0]);
        snap.y.push_back(v[1]);
        snap.z.push_back(v[2]);
        snap.vx.push_back(v[3]);
        snap.vy.push_back(v[4]);
        snap.vz.push_back(v[5]);
        snap.m.push_back(v[6]);
        snap.h.push_back(v[7]);
        snap.rho.push_back(v[8]);
    }
    if (snap.size() == 0)
        throw std::runtime_error(
            "SPH snapshot file (plain_ascii format) contained no rows -- expected 9 columns "
            "(x y z vx vy vz m h rho) per line: " + path);
    return snap;
}

// `splash to ascii`'s preamble is a fixed shape (signature line, then
// time/npart/units blocks, each one `#` comment line followed by exactly
// one uncommented values line, then a final `#` column-header line) but we
// deliberately don't hardcode its length or the exact number of blocks --
// npart in particular carries one column per particle type SPLASH knows
// about, which has grown over the years and could again. Instead: scan
// forward, skipping `#` comment lines; whenever an uncommented line is
// found, peek at the line after it -- if THAT line is also uncommented
// (and non-blank), the first line is the start of the per-particle data
// block (every subsequent row is uncommented too); otherwise the first
// line is a single metadata values line (time/npart/units) sandwiched
// between comments, so skip it and keep scanning. This finds the data
// block regardless of how many metadata blocks precede it.
//
// Returns the 0-based line number the data block starts at.
int find_splash_ascii_data_start(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open SPH snapshot file: " + path);

    std::string cur;
    if (!std::getline(f, cur))
        throw std::runtime_error("splash-ascii snapshot file is empty: " + path);

    int line_no = 0;
    while (true) {
        if (is_comment_line(cur) || is_blank_line(cur)) {
            ++line_no;
            if (!std::getline(f, cur))
                throw std::runtime_error(
                    "splash-ascii snapshot: reached end of file while still in the header "
                    "preamble, no particle data found: " + path);
            continue;
        }
        // `cur` is an uncommented line -- peek ahead.
        std::string next;
        bool have_next = static_cast<bool>(std::getline(f, next));
        if (have_next && !is_comment_line(next) && !is_blank_line(next)) {
            return line_no; // `cur` is the first particle-data row
        }
        // `cur` was a single metadata values line (time/npart/units); skip it.
        ++line_no;
        if (!have_next)
            throw std::runtime_error(
                "splash-ascii snapshot: reached end of file while still in the header "
                "preamble, no particle data found: " + path);
        cur = next;
    }
}

int count_whitespace_tokens(const std::string& s) {
    std::istringstream iss(s);
    int n = 0;
    std::string tok;
    while (iss >> tok) ++n;
    return n;
}

// Strips a bracketed unit suffix (e.g. " [pc]", " [M_{Sun}]", " [km/s]") --
// there's at most one per label in real `splash to ascii` output -- and
// collapses whitespace/underscores/backslashes, so "particle mass
// [M_{Sun}]", "v_x [km/s]", and "\alpha [1/s]" all normalize down to
// "particlemass", "vx", and "alpha" respectively. This lets us match a
// label to a semantic column identity independent of exactly how SPLASH
// chose to format it.
std::string normalize_splash_label(const std::string& raw) {
    std::string s;
    bool in_brackets = false;
    for (char c : raw) {
        if (c == '[') { in_brackets = true; continue; }
        if (c == ']') { in_brackets = false; continue; }
        if (in_brackets) continue;
        if (c == '\\' || c == '_' || c == ' ' || c == '\t') continue;
        s += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

// Extracts the unit string inside the LAST "[...]" bracket pair of a raw
// (un-normalized) splash_ascii column label -- e.g. "particle mass
// [M_{Sun}]" -> "M_{Sun}", "density [g/cm^3]" -> "g/cm^3". Returns "" if
// the label has no bracket at all.
std::string extract_label_unit(const std::string& raw) {
    size_t close = raw.rfind(']');
    if (close == std::string::npos) return "";
    size_t open = raw.rfind('[', close);
    if (open == std::string::npos) return "";
    return raw.substr(open + 1, close - open - 1);
}

std::string normalize_unit(const std::string& raw) {
    std::string s;
    for (char c : raw) {
        if (c == '{' || c == '}' || c == '_' || c == ' ' || c == '\t') continue;
        s += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

bool is_msun_unit(const std::string& raw) { return normalize_unit(raw) == "msun"; }
bool is_pc_unit(const std::string& raw) { return normalize_unit(raw) == "pc"; }
bool is_gcm3_unit(const std::string& raw) {
    std::string u = normalize_unit(raw);
    return u == "g/cm^3" || u == "g/cm3" || u == "gcm^-3" || u == "gcm-3";
}

// g/cm^3 per Msun/pc^3, using the same constants as this project's own
// SPLASH-comparison scripts (compare_splash.py, investigate_discrepancy.py,
// verify_h_hypothesis.py) -- kept in one place so the h<->density
// derivation below and those scripts can never quietly drift apart.
double msun_pc3_to_g_cm3() {
    constexpr double msun_g = 1.9891e33;
    constexpr double pc_cm = 3.08567758e18;
    return msun_g / (pc_cm * pc_cm * pc_cm);
}

// Only the columns KIWI actually needs -- or needs to recognize in
// order to treat specially, like itype -- are listed here. Any other
// column present in a given file (pressure, alpha, div_v, poten, dt, or a
// physics-module column this codebase has never heard of -- Bx/By/Bz for
// an MHD run, dust fraction, radiation quantities, etc.) is simply
// ignored rather than rejected: `splash to ascii` writes whatever columns
// the source Phantom dump actually has, in whatever order, and that set
// genuinely varies run to run (see write_sphdump() in SPLASH's own
// write_sphdata.f90 -- ncolumns and the column order are both dynamic,
// and even the presence of the itype column depends on whether the dump
// has more than one particle type at all).
const std::unordered_map<std::string, std::string>& splash_ascii_column_aliases() {
    static const std::unordered_map<std::string, std::string> aliases = {
        {"x", "x"}, {"y", "y"}, {"z", "z"},
        {"particlemass", "m"}, {"mass", "m"}, {"m", "m"},
        {"h", "h"}, {"smoothinglength", "h"},
        {"density", "rho"}, {"rho", "rho"},
        {"vx", "vx"}, {"vy", "vy"}, {"vz", "vz"},
        {"itype", "itype"},
    };
    return aliases;
}

// Splits a splash_ascii column-header comment line ("# x [pc]  y [pc] ...
// itype") into one label string per data column, in file order.
//
// This can't be done with a plain whitespace split: `splash to ascii`
// writes each label left-justified in a fixed-width field (24 or 16
// characters -- 23 or 15 content columns plus a 1-space separator,
// depending on whether SPLASH itself was built in double or single
// precision; see write_sphdump(), src/write_sphdata.f90) with no other
// delimiter, so a multi-word label like "particle mass [M_{Sun}]" or "div
// v [erg]" can't be told apart from its neighbors by whitespace alone.
// The one exception is a trailing `itype` column: SPLASH appends it with
// no field-width padding at all (it's passed to the write statement as a
// separate scalar argument, not part of the per-column label array), and
// only when the dump has more than one particle type -- a pure-gas dump
// with no sinks/dust/stars gets no itype column in the file at all.
//
// `n_data_cols` is the column count already known unambiguously from a
// whitespace split of an actual data row (numbers never contain embedded
// spaces, so that split is never ambiguous) -- this function trusts that
// count rather than trying to independently infer it from the header.
std::vector<std::string> split_splash_ascii_labels(const std::string& header_line,
                                                     int n_data_cols,
                                                     const std::string& path) {
    if (header_line.size() < 2 || header_line[0] != '#')
        throw std::runtime_error(
            "splash-ascii snapshot: expected the line immediately before the particle data "
            "to be the '#'-prefixed column-header line, but it wasn't: " + path);
    std::string body = header_line.substr(1);
    while (!body.empty() && (body.back() == '\r' || body.back() == '\n')) body.pop_back();
    if (!body.empty() && body.front() == ' ') body = body.substr(1); // the format's `1x` after `#`

    // SPLASH always writes the literal lowercase token `itype` with no
    // padding after it, iff there's more than one particle type -- and
    // nothing else in a real column label ends that way, so this is safe
    // to detect by suffix regardless of any trailing-whitespace stripping
    // a text editor or version-control filter might have done elsewhere
    // on the line (trailing-whitespace stripping only removes whitespace,
    // never the last non-whitespace token).
    bool has_itype_tail = body.size() >= 5 && body.compare(body.size() - 5, 5, "itype") == 0;

    int n_fixed = has_itype_tail ? n_data_cols - 1 : n_data_cols;
    if (n_fixed <= 0)
        throw std::runtime_error(
            "splash-ascii snapshot: column-header line doesn't have enough fields for the " +
            std::to_string(n_data_cols) + " column(s) seen in the particle data: " + path);

    std::string fixed_part = has_itype_tail ? body.substr(0, body.size() - 5) : body;

    // Field width is one of exactly two values, fixed by SPLASH's own
    // source (23 or 15 content characters + a 1-space separator) -- not
    // something to re-derive by dividing this particular line's length,
    // which would silently misalign every field if the *last* field's
    // trailing padding happened to be stripped (a real risk when there's
    // no itype tail holding the line length steady). Try both candidates;
    // when there IS an itype tail, `fixed_part` is never truncated (itype
    // comes right after it, not end-of-line), so the division must land
    // exactly; when there isn't, allow the final field to be short.
    int field_w = -1;
    for (int candidate : {24, 16}) {
        if (has_itype_tail) {
            if (fixed_part.size() == static_cast<size_t>(n_fixed) * candidate) { field_w = candidate; break; }
        } else {
            size_t min_len = static_cast<size_t>(n_fixed - 1) * candidate + 1;
            size_t max_len = static_cast<size_t>(n_fixed) * candidate;
            if (fixed_part.size() >= min_len && fixed_part.size() <= max_len) { field_w = candidate; break; }
        }
    }
    if (field_w < 0)
        throw std::runtime_error(
            "splash-ascii snapshot: could not reliably split the column-header line into " +
            std::to_string(n_fixed) + " field(s) -- its length doesn't match either of "
            "SPLASH's known label field widths (23 or 15 characters): " + path);

    std::vector<std::string> labels;
    labels.reserve(static_cast<size_t>(n_data_cols));
    for (int i = 0; i < n_fixed; ++i) {
        std::string field = (!has_itype_tail && i == n_fixed - 1)
            ? fixed_part.substr(static_cast<size_t>(i) * field_w) // rest of the string --
                                                                    // tolerant of stripped
                                                                    // trailing padding
            : fixed_part.substr(static_cast<size_t>(i) * field_w,
                                 static_cast<size_t>(field_w));
        size_t a = field.find_first_not_of(" \t");
        size_t b = field.find_last_not_of(" \t");
        labels.push_back(a == std::string::npos ? "" : field.substr(a, b - a + 1));
    }
    if (has_itype_tail) labels.push_back("itype");
    return labels;
}

SPHSnapshot load_splash_ascii(const std::string& path, double hfact) {
    int data_start = find_splash_ascii_data_start(path);
    if (data_start == 0)
        throw std::runtime_error(
            "splash-ascii snapshot: no header preamble found before the particle data -- "
            "expected at least the column-header comment line: " + path);

    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open SPH snapshot file: " + path);
    std::string line, header_line;
    for (int i = 0; i < data_start; ++i) {
        if (!std::getline(f, line))
            throw std::runtime_error("splash-ascii snapshot: file changed size while reading "
                                      "the header preamble: " + path);
        header_line = line; // after the loop: the line right before the data -- the
                             // column-header comment row
    }
    std::string first_data_line;
    if (!std::getline(f, first_data_line))
        throw std::runtime_error("splash-ascii snapshot: no particle data found: " + path);

    int n_data_cols = count_whitespace_tokens(first_data_line);
    if (n_data_cols < 8) // sanity floor: fewer than x,y,z,m,vx,vy,vz + (h or rho) can't be real
        throw std::runtime_error(
            "splash-ascii snapshot: first particle-data row has only " +
            std::to_string(n_data_cols) + " column(s), too few to be a real splash_ascii "
            "file: " + path);

    std::vector<std::string> labels = split_splash_ascii_labels(header_line, n_data_cols, path);

    std::unordered_map<std::string, int> col_index; // semantic name -> column index
    const auto& aliases = splash_ascii_column_aliases();
    for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
        auto it = aliases.find(normalize_splash_label(labels[i]));
        if (it == aliases.end()) continue; // an unrecognized/unneeded column -- ignore it
        auto ins = col_index.emplace(it->second, i);
        if (!ins.second)
            throw std::runtime_error(
                "splash-ascii snapshot: the column-header line appears to name '" + it->second +
                "' more than once (columns " + std::to_string(ins.first->second) + " and " +
                std::to_string(i) + ") -- the header likely failed to parse correctly: " + path);
    }

    for (const std::string& req : {std::string("x"), std::string("y"), std::string("z"),
                                    std::string("m"), std::string("vx"), std::string("vy"),
                                    std::string("vz")}) {
        if (!col_index.count(req))
            throw std::runtime_error(
                "splash-ascii snapshot: could not find a '" + req + "' column in the "
                "column-header line: " + path);
    }
    bool have_h = col_index.count("h") != 0;
    bool have_rho = col_index.count("rho") != 0;
    if (!have_h && !have_rho)
        throw std::runtime_error(
            "splash-ascii snapshot: could not find either an 'h' (smoothing length) or a "
            "'density' column in the column-header line -- need at least one, to derive the "
            "other from it: " + path);
    bool have_itype = col_index.count("itype") != 0;

    if (!have_itype) {
        std::cerr << "[kiwi] splash-ascii snapshot: no 'itype' column in the header -- this file "
                     "has only one particle type (SPLASH omits itype entirely in that case); "
                     "treating every row as gas: " << path << "\n";
    }
    if (have_h != have_rho) {
        // h = hfact*(m/rho)^(1/3), hfact = 1.2 for the cubic spline -- Price et al. 2018 (the Phantom paper), eq. 6.
        // hfact is the caller's (Config::kernel_hfact(): per-kernel Phantom default, or kernel.hfact).
        // Note this is the *inverse* of (m/rho)^(-1/3): that exponent sign would make h grow
        // with density, which is backwards (denser gas is more highly resolved -> smaller h).
        //
        // That formula only holds for m, h, and rho expressed in a mutually consistent unit
        // system -- but `splash to ascii` writes m in M_Sun, h in pc, and density in g/cm^3
        // (three DIFFERENT unit systems), so deriving one from the other requires converting
        // through Msun/pc^3 first. Rather than assume every file uses exactly those units,
        // check the column-header line's own unit brackets and refuse (loudly) if they say
        // anything else, instead of silently producing numbers in the wrong units.
        std::string m_unit = extract_label_unit(labels[static_cast<size_t>(col_index.at("m"))]);
        std::string present_unit = extract_label_unit(
            labels[static_cast<size_t>(col_index.at(have_h ? "h" : "rho"))]);
        bool units_ok = is_msun_unit(m_unit) &&
            (have_h ? is_pc_unit(present_unit) : is_gcm3_unit(present_unit));
        if (!units_ok)
            throw std::runtime_error(
                "splash-ascii snapshot: only a '" + std::string(have_h ? "h" : "density") +
                "' column was found, and this tool only knows how to derive the other from it "
                "when mass is in M_Sun, h is in pc, and density is in g/cm^3 (the h = "
                "hfact*(m/rho)^(1/3) relation needs m, h, and rho in a mutually consistent unit "
                "system) -- got mass unit '" + m_unit + "' and " +
                (have_h ? "h" : "density") + " unit '" + present_unit + "' instead. Extend "
                "this derivation before using it on a file in different units: " + path);
        std::cerr << "[kiwi] splash-ascii snapshot: only a '" << (have_h ? "h" : "density")
                  << "' column was found in the header -- deriving "
                  << (have_h ? "density from h" : "h from density")
                  << " via h = hfact*(m/rho)^(1/3) with hfact = " << hfact << " (Price et al. 2018, eq. 6), converting "
                     "through Msun/pc^3 since density is in g/cm^3 but mass and h aren't: "
                  << path << "\n";
    }

    int idx_x = col_index.at("x"), idx_y = col_index.at("y"), idx_z = col_index.at("z");
    int idx_m = col_index.at("m");
    int idx_vx = col_index.at("vx"), idx_vy = col_index.at("vy"), idx_vz = col_index.at("vz");
    int idx_h = have_h ? col_index.at("h") : -1;
    int idx_rho = have_rho ? col_index.at("rho") : -1;
    int idx_itype = have_itype ? col_index.at("itype") : -1;

    SPHSnapshot snap;
    long long n_rows_seen = 0, n_dropped_non_gas = 0;
    int row = 0;
    std::vector<double> vals(static_cast<size_t>(n_data_cols));

    auto process_row = [&](const std::string& row_line) {
        if (is_blank_line(row_line)) return;
        if (is_comment_line(row_line)) return; // tolerate a trailing/footer comment, if any

        std::istringstream iss(row_line);
        for (int i = 0; i < n_data_cols; ++i) {
            if (!(iss >> vals[static_cast<size_t>(i)]))
                throw std::runtime_error(
                    "splash-ascii snapshot: could not parse a particle row (expected " +
                    std::to_string(n_data_cols) + " columns, matching the column-header "
                    "line): " + path + ": " + row_line);
        }

        ++n_rows_seen;
        long itype = have_itype ? std::lround(vals[static_cast<size_t>(idx_itype)]) : 1L;
        if (itype != 1) {
            ++n_dropped_non_gas;
            return;
        }

        double m = vals[static_cast<size_t>(idx_m)];
        double h, rho;
        if (have_h && have_rho) {
            h = vals[static_cast<size_t>(idx_h)];
            rho = vals[static_cast<size_t>(idx_rho)];
        } else if (have_h) {
            // m is in Msun and h is in pc (checked above), but rho must come out in g/cm^3 to
            // match every other splash_ascii density column -- so derive in Msun/pc^3 first
            // (the unit system the hfact relation actually needs) and convert at the end.
            h = vals[static_cast<size_t>(idx_h)];
            double ratio = h / hfact;
            double rho_msun_pc3 = (ratio > 0.0) ? m / (ratio * ratio * ratio) : 0.0;
            rho = rho_msun_pc3 * msun_pc3_to_g_cm3();
        } else {
            // Inverse: density column is in g/cm^3, so convert to Msun/pc^3 before applying
            // h = hfact*(m/rho)^(1/3) (which needs m and rho in a consistent unit system), so
            // the resulting h comes out in pc, matching every other splash_ascii h column.
            rho = vals[static_cast<size_t>(idx_rho)];
            double rho_msun_pc3 = rho / msun_pc3_to_g_cm3();
            h = (m > 0.0 && rho_msun_pc3 > 0.0) ? hfact * std::cbrt(m / rho_msun_pc3) : 0.0;
        }

        ++row;
        snap.id.push_back(row);
        snap.x.push_back(vals[static_cast<size_t>(idx_x)]);
        snap.y.push_back(vals[static_cast<size_t>(idx_y)]);
        snap.z.push_back(vals[static_cast<size_t>(idx_z)]);
        snap.vx.push_back(vals[static_cast<size_t>(idx_vx)]);
        snap.vy.push_back(vals[static_cast<size_t>(idx_vy)]);
        snap.vz.push_back(vals[static_cast<size_t>(idx_vz)]);
        snap.m.push_back(m);
        snap.h.push_back(h);
        snap.rho.push_back(rho);
    };

    process_row(first_data_line);
    while (std::getline(f, line)) process_row(line);

    if (n_dropped_non_gas > 0) {
        std::cerr << "[kiwi] splash-ascii snapshot: dropped " << n_dropped_non_gas << " of "
                  << n_rows_seen << " particle row(s) with itype != 1 (non-gas -- sink/dust/"
                     "star/etc.); only gas particles are used for the grid: " << path << "\n";
    }
    if (snap.size() == 0)
        throw std::runtime_error(
            "splash-ascii snapshot file contained no gas (itype == 1) particle rows: " + path);
    return snap;
}

} // namespace

std::string detect_snapshot_format(const std::string& path) {
    if (is_phantom_dump(path)) return "phantom_binary";
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open SPH snapshot file: " + path);
    std::string first_line;
    std::getline(f, first_line);
    // splash-to-ascii files always open with a `#` comment naming SPLASH as
    // the tool that wrote them (see the module header comment for the full
    // preamble shape); the plain_ascii format has no header at all --
    // its first line is already the first particle's data.
    if (!first_line.empty() && first_line[0] == '#' && first_line.find("SPLASH") != std::string::npos) {
        return "splash_ascii";
    }
    return "plain_ascii";
}

SPHSnapshot load_sph_snapshot(const std::string& path, const std::string& format, double hfact,
                              const PhantomUnits& units, bool allow_missing_velocity) {
    std::string resolved = format;
    if (resolved == "auto") resolved = detect_snapshot_format(path);

    if (resolved == "plain_ascii") return load_plain_ascii(path);
    if (resolved == "splash_ascii") return load_splash_ascii(path, hfact);
    if (resolved == "phantom_binary") return load_phantom_dump(path, units, allow_missing_velocity);
    throw std::runtime_error(
        "unknown SPH snapshot format '" + format + "' (expected auto, plain_ascii, "
        "splash_ascii, or phantom_binary)");
}

void write_plain_ascii(const SPHSnapshot& snap, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("could not open output file for writing: " + path);
    f.precision(17);
    for (size_t i = 0; i < snap.size(); ++i) {
        f << snap.x[i] << ' ' << snap.y[i] << ' ' << snap.z[i] << ' ' << snap.vx[i] << ' ' << snap.vy[i] << ' '
          << snap.vz[i] << ' ' << snap.m[i] << ' ' << snap.h[i] << ' ' << snap.rho[i] << '\n';
    }
    if (!f) throw std::runtime_error("error while writing: " + path);
}

PointSet extract_points(const SPHSnapshot& snap) {
    PointSet points;
    points.id = snap.id;
    points.x = snap.x;
    points.y = snap.y;
    points.z = snap.z;
    return points;
}

} // namespace kiwi
