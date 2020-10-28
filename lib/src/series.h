#ifndef SERIES_H
#define SERIES_H

#include "gf4.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <functional>

namespace cppcodes{
class Series{
    public:
        std::vector<gf4> coeffs;
        size_t zero_shift;
        Series(): coeffs(), zero_shift(0){ coeffs.push_back(gf4()); };
        Series(const Series& o): coeffs(o.coeffs), zero_shift(o.zero_shift){};
        Series(Series&& o): coeffs(std::move(o.coeffs)), zero_shift(o.zero_shift){};
        Series(size_t n, size_t shift): coeffs(n), zero_shift(shift) 
        { if (zero_shift < 0) throw std::invalid_argument("Negative shift"); };
        Series(std::vector<gf4> coeffs_, size_t shift = 0): coeffs(coeffs_), zero_shift(shift) 
        { if (zero_shift < 0) throw std::invalid_argument("Negative shift"); };
        Series(std::string coeffs_, size_t shift = 0) {
            zero_shift = shift;
            coeffs.reserve(coeffs_.size());
            for (size_t i = 0; i < coeffs_.size(); ++i)
                coeffs.push_back(gf4(coeffs_[i]));
        }

        Series& operator=(const Series& o){
            if (this != &o){
                zero_shift = o.zero_shift;
                coeffs.assign(o.coeffs.begin(), o.coeffs.end());
            }
            return *this;
        };

        gf4& operator[](size_t i){
            if (zero_shift + i < 0 || zero_shift +i >= coeffs.size())
                throw std::out_of_range("Index is out of bound");
            return coeffs[zero_shift + i];
        };

        const gf4& operator[](size_t i) const{
            if (zero_shift + i < 0 || zero_shift +i >= coeffs.size())
                throw std::out_of_range("Index is out of bound");
            return coeffs[zero_shift + i];      
        };

        gf4 at(size_t i) const{
            if (zero_shift + i < 0 || zero_shift + i >= coeffs.size())
                return gf4();
            return coeffs[zero_shift + i];
        };

        int min_power() const {
            return -zero_shift;
        };

        int max_power() const {
            return coeffs.size() - zero_shift - 1;
        };

        Series& strip(){
            if (coeffs.size() == 1)
                return *this;
            
            size_t i = 0;
            while(i < coeffs.size() && coeffs[i] == 0 && zero_shift > 0)
            {
                ++i;
                --zero_shift;
            }
            if (i > 0){
                for (size_t j = 0; j < coeffs.size() - i; ++j)
                    coeffs[j] = coeffs[j + i];
                while (i-- > 0)
                    coeffs.pop_back();
            }
            while (coeffs.size() > zero_shift + 1 && coeffs.back() == 0)
                coeffs.pop_back();

            return *this;
        };

        Series operator+(const Series& b) const {
            size_t c_shift = 0;
            std::vector<gf4> ks;
            int right = std::max(max_power(), b.max_power());
            bool first = true;
            for (int i=std::min(b.min_power(), min_power()); i <= right; ++i){
                gf4 s = at(i) + b.at(i);
                if (!first || !(s == 0) || i == 0)
                {
                    ks.push_back(s);
                    if (first)
                    {
                        c_shift = -i;
                        first = false;
                    }
                }
            }
            while (ks.size() > (c_shift + 1) && ks.back() == 0)
            {
                ks.pop_back();
            }
            return Series(ks, c_shift);
        };

        Series operator*(const Series& b) const {
            size_t c_shift = zero_shift + b.zero_shift;
            int len = coeffs.size() + b.coeffs.size() - 1;
            std::vector<gf4> ks(len);
            for (size_t i = 0; i < coeffs.size(); ++i)
                for (size_t j = 0; j < b.coeffs.size(); ++j)
                    ks[i + j] = ks[i + j] + coeffs[i] * b.coeffs[j];
            Series s(ks, c_shift);
            s.strip();
            return s;
        };

        bool operator<(const Series& o){
            if (this == &o){
                return false;
            }
            if (zero_shift < o.zero_shift){
                return true;
            }
            for (size_t i = 0; i < coeffs.size() && i < o.coeffs.size(); ++i)
                if (coeffs[i] != o.coeffs[i]){
                    return coeffs[i].value < o.coeffs[i].value;
                }
            return coeffs.size() < o.coeffs.size();
        }

        Series inverse() const {
            Series s(*this);
            std::reverse(s.coeffs.begin(), s.coeffs.end());
            s.zero_shift = s.coeffs.size() - 1 - s.zero_shift;
            s.strip();
            return s;
        };

        Series conj() const{
            std::vector<gf4> ks;
            ks.reserve(coeffs.size());
            for (size_t i = 0; i < coeffs.size(); ++i)
                ks.push_back(coeffs[i].conj());
            return Series(ks, zero_shift);
        };

        gf4 dotProduct(const std::vector<gf4>& v){
            gf4 s;
            size_t j = 0;
            for (int i = v.size() - 1; i >= 0; --i){
                s = s + (v[i] * at(j++));
            }
            return s;
        }

        std::string toString() const {
            std::string repr = "";
            for (size_t i = 0; i < coeffs.size(); ++i)
                if (i == zero_shift && zero_shift != 0){
                    repr.append("(");
                    repr.append(coeffs[i].toString());
                    repr.append(")");
                } else
                    repr.append(coeffs[i].toString());
            
            return repr;
        }

        bool operator==(const Series& o) const {
            if (this == &o)
                return true;
            
            if (zero_shift != o.zero_shift || coeffs.size() != o.coeffs.size())
                return false;
            for (size_t i = 0; i < coeffs.size(); ++i)
                if (coeffs[i] != o.coeffs[i])
                    return false;
            
            return true;
        }
};

struct SeriesHasher{
    std::size_t operator()(const Series& s) const {
        int v = 0;
        for (auto it = s.coeffs.begin(); it != s.coeffs.end(); ++it)
            v = v * GF4_SIZE + it->value;
        return std::hash<int>{}((int)s.zero_shift) ^ std::hash<int>{}(v);
    }
};
}

#endif