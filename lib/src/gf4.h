#ifndef GF4_H
#define GF4_H

#include <string>
#include <stdexcept>

const short PROD[][4] = {
    {0, 0, 0, 0},
    {0, 1, 2, 3},
    {0, 2, 3, 1},
    {0, 3, 1, 2},
};

const short SUM[][4] = {
    {0, 1, 2, 3},
    {1, 0, 3, 2},
    {2, 3, 0, 1},
    {3, 2, 1, 0},
};

const std::string REPR[] = {"0", "1", "u", "v"};

const size_t GF4_SIZE = 4;

namespace cppcodes{
class gf4{
    public: 
        short value;
        gf4(): value(0) {};
        gf4(char v) {
            if (v >= 0 && v <= 4)
                value = v;
            else{
                for (size_t i = 0; i < sizeof(REPR); ++i)
                if (v == REPR[i][0]){
                    value = i;
                    return;
                }
                throw std::invalid_argument("unrecognized value");
            }
        };
        gf4(std::string v) {
            for (size_t i = 0; i < sizeof(REPR); ++i)
                if (v == REPR[i]){
                    value = i;
                    return;
                }
            throw std::invalid_argument("unrecognized value");
        };

        short trace() const { return value < 2 ? 0 : 1; };
        gf4 conj() const { return value < 2 ? gf4(value) : gf4(5 - value); }
        std::string toString() const { return REPR[value]; };

        gf4& operator++(){ value = (value + 1) % GF4_SIZE; return *this; }
        gf4 operator*(const gf4& b) const { return gf4(PROD[value][b.value]);} 
        gf4 operator+(const gf4& b) const {return gf4(SUM[value][b.value]); };
        bool operator==(const gf4& b) const { return value == b.value; };
        bool operator!=(const gf4& b) const { return value != b.value; }
        bool operator==(const char b) const { return value == b; }
        bool operator!=(const char b) const { return value != b; }
        bool operator==(const std::string& b) const {
            for (size_t i = 0; i < sizeof(REPR); ++i)
                if (b == REPR[i])
                    return value == (short)i;
            return false;
        }
        bool operator!=(const std::string& b) const { return !(*this == b); }
};

}

#endif
