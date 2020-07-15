#ifndef CODES_H
#define CODES_H

#include "series.h"
#include "gf4.h"
#include <vector>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <queue>
#include <climits>
#include <memory>

#include <iostream>

typedef long long int xlong;
const long long int INFTY = LLONG_MAX;

namespace cppcodes{

class Combinations{
    size_t n;
    size_t k;
    std::shared_ptr<std::vector<std::vector<size_t>>> results; 
    std::vector<size_t> current;
    
    void generate(size_t i, size_t j);

    public:
        Combinations(size_t n_, size_t k_): n(n_), k(k_) {}
        std::shared_ptr<std::vector<std::vector<size_t>>> get();    
};

class Code{
    public:
        std::vector<Series> generators;
        size_t n;
        size_t k;
        Code(std::vector<Series> gens_, size_t n_, size_t k_ = 1): generators(gens_), n(n_), k(k_) {};
        Code(size_t n_, size_t k_ = 1): generators(), n(n_), k(k_) {};
        
        Series& operator[](size_t i);
        const Series& operator[](size_t i) const;
        void add(Series& g);
        void remove(size_t i);
        bool validate();
        size_t maxSize(size_t i);
        bool isSelfOrthogonal();
        bool isOrthogonal(Code& other);
        std::string toString();
        xlong minDistance();
};

class CodeGenerator {
    size_t n;
    size_t k;
    size_t nu;
    size_t current_degree_split;
    std::vector<gf4> current_code;
    std::shared_ptr<std::vector<std::vector<size_t>>> degrees;

    public:
    CodeGenerator(size_t n_, size_t k_, size_t nu_);
    Code next();
};

}

#endif