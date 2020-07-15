#ifndef SEARCH_H
#define SEARCH_H

#include <stdexcept>
#include <memory>
#include "codes.h"
#include "series.h"
#include "gf4.h"

namespace cppcodes{
class SearchSelfOrthogonal{
private:
    bool warm;
    size_t n;
    size_t degree;
    size_t k;
    std::unordered_multimap<Series, std::shared_ptr<std::vector<Series>>, SeriesHasher> Rgg;
    std::vector<Code> codes;
public:
    SearchSelfOrthogonal(size_t n_, size_t degree_, size_t k_ = 1)
    : warm(false)
    , n(n_)
    , degree(degree_)
    , k(k_)
    , Rgg()
    , codes() {
        if (k != 1)
            throw new std::logic_error("Not implemented.");
    };

    size_t getN(){
        return n;
    }

    size_t getK(){
        return k;
    }

    size_t getDegree(){
        return degree;
    }

    void initialize(){
        if (warm) return;
        std::vector<gf4> s(degree + 1);
        s[0] = gf4(1);
        for(;;){
            Series series(s);
            auto v = (series.inverse().conj() * series).strip();
            auto search = Rgg.find(v);
            if (search == Rgg.end()){
                Rgg.insert({v, std::shared_ptr<std::vector<Series>>(new std::vector<Series>())});
            }
            search = Rgg.find(v);
            search->second->push_back(series);
            // find next series
            bool carry = true;
            for (int i = degree; i >= 0; --i)
                if (carry){
                    carry = ++s[i] == 0; 
                }
            if (s[0] != 1)
                break;
        }
        warm = true;
    }

    void append(const std::vector<Series>& rgg, std::vector<Series>& code, size_t i){
        if (i == n) {
            //codes.push_back(Code(code, n, k));
        } else {
            auto v = Rgg.find(rgg[i])->second;
            for (auto it = v->begin(); it != v->end(); ++it){
                code.push_back(*it);
                append(rgg, code, i + 1);
                code.pop_back();
            }
        }
    }

    void generate(std::vector<Series>& v, const Series& s, size_t i){
        if (i == n - 1){
            auto search = Rgg.find(s);
            if (search != Rgg.end()){
                v.push_back(s);
                std::vector<Series> code;
                code.reserve(n);
                append(v, code, 0);
                v.pop_back();
            }
        } else {
            for (auto it = Rgg.begin(); it != Rgg.end(); ++it){
                v.push_back(it->first);
                generate(v, s + it->first, i + 1);
                v.pop_back();
            }
        }
    }

    std::vector<Code> find(){
        initialize();
        std::vector<Series> v;
        v.reserve(n);
        Series s;
        generate(v, s, 0);
        return codes;
    }
};
}

#endif