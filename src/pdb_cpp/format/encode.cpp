#include "encode.h"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

const std::string kDigitsUpper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const std::string kDigitsLower = "0123456789abcdefghijklmnopqrstuvwxyz";

long long pow_int(long long base, int exp) {
    long long result = 1;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

bool is_all_spaces(const std::string& s) {
    for (char c : s) {
        if (c != ' ') {
            return false;
        }
    }
    return true;
}

long long parse_decimal_like_python(const std::string& s) {
    size_t idx = 0;
    long long value = std::stoll(s, &idx, 10);
    for (size_t i = idx; i < s.size(); ++i) {
        if (!std::isspace(static_cast<unsigned char>(s[i]))) {
            throw std::invalid_argument("invalid number literal.");
        }
    }
    return value;
}

int decode_pure_upper(const std::string& s) {
    int result = 0;
    for (char c : s) {
        result *= 36;
        if (c >= '0' && c <= '9') {
            result += c - '0';
        } else if (c >= 'A' && c <= 'Z') {
            result += c - 'A' + 10;
        } else {
            throw std::invalid_argument("invalid number literal.");
        }
    }
    return result;
}

int decode_pure_lower(const std::string& s) {
    int result = 0;
    for (char c : s) {
        result *= 36;
        if (c >= '0' && c <= '9') {
            result += c - '0';
        } else if (c >= 'a' && c <= 'z') {
            result += c - 'a' + 10;
        } else {
            throw std::invalid_argument("invalid number literal.");
        }
    }
    return result;
}

int decode_basic_hex_cap(const std::string& s) {
    int result = 0;
    for (char c : s) {
        result *= 16;
        if (c == ' ') {
            result += 0;
        } else if (c >= '0' && c <= '9') {
            result += c - '0';
        } else if (c >= 'A' && c <= 'Z') {
            result += c - 'A' + 10;
        } else {
            throw std::invalid_argument("invalid number literal.");
        }
    }
    return result;
}

std::string encode_pure(const std::string& digits, long long value) {
    if (value < 0) {
        throw std::invalid_argument("value out of range.");
    }
    if (value == 0) {
        return digits.substr(0, 1);
    }
    const long long base = static_cast<long long>(digits.size());
    std::string result;
    while (value != 0) {
        long long rest = value / base;
        result.push_back(digits[static_cast<size_t>(value - rest * base)]);
        value = rest;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

} // namespace

int hy36decode(int width, const std::string& s) {
    if (static_cast<int>(s.size()) != width) {
        throw std::invalid_argument("invalid number literal.");
    }

    char f = s[0];
    if (f == '-' || f == ' ' || std::isdigit(static_cast<unsigned char>(f))) {
        try {
            return static_cast<int>(parse_decimal_like_python(s));
        } catch (const std::exception&) {
            try {
                return decode_basic_hex_cap(s);
            } catch (const std::exception&) {
                // fall through
            }
        }
        if (is_all_spaces(s)) {
            return 0;
        }
    } else if (f >= 'A' && f <= 'Z') {
        long long decoded = decode_pure_upper(s);
        long long offset = 10 * pow_int(36, width - 1);
        long long result = decoded - offset + pow_int(10, width);
        return static_cast<int>(result);
    } else if (f >= 'a' && f <= 'z') {
        long long decoded = decode_pure_lower(s);
        long long offset = 16 * pow_int(36, width - 1);
        long long result = decoded + offset + pow_int(10, width);
        return static_cast<int>(result);
    }

    throw std::invalid_argument("invalid number literal.");
}

std::string hy36encode(int width, int value) {
    long long i = value;
    long long min_value = 1 - pow_int(10, width - 1);
    if (i < min_value) {
        throw std::invalid_argument("value out of range.");
    }

    long long pow10 = pow_int(10, width);
    if (i < pow10) {
        std::ostringstream oss;
        oss << std::setw(width) << i;
        return oss.str();
    }

    i -= pow10;
    long long block = 26 * pow_int(36, width - 1);
    if (i < block) {
        i += 10 * pow_int(36, width - 1);
        return encode_pure(kDigitsUpper, i);
    }

    i -= block;
    if (i < block) {
        i += 10 * pow_int(36, width - 1);
        return encode_pure(kDigitsLower, i);
    }

    throw std::invalid_argument("value out of range.");
}
