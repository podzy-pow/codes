#include "codes.h"

#define LOG(msg) \
    std::cout << __FILE__ << "(" << __LINE__ << "): " << msg << std::endl 

using namespace cppcodes;

void Combinations::generate(size_t i, size_t j){
    if (j > i)
        return;
    if (i == 0){
        results->push_back(current);
        return;
    }
    if (j == 0)
        return;
    for (size_t k = 1; k <= (i - j + 1); ++k){
        current.push_back(k);
        generate(i - k, j - 1);
        current.pop_back();
    }
}

std::shared_ptr<std::vector<std::vector<size_t>>> Combinations::get() {
    results = std::make_shared<std::vector<std::vector<size_t>>>();
    current.clear();
    if (k <= n)
        generate(n, k);
    return results;
}

Series& Code::operator[](size_t i) {
    if (i < 0 || i >= generators.size())
        throw std::out_of_range("Index is out of bound");
    return generators[i];
}

const Series& Code::operator[](size_t i) const {
    if (i < 0 || i >= generators.size())
        throw std::out_of_range("Index is out of bound");
    return generators[i];
}

void Code::add(Series& g){
    Series s(g);
    generators.push_back(s);
}

void Code::remove(size_t i){
    if (i < 0 || i >= generators.size())
        throw std::out_of_range("Index is out of bound");
    generators.erase(generators.begin() + i);
}

bool Code::validate() {
    bool ok = true;
    for (auto gen: generators){
        ok = ok && gen.zero_shift == 0 && gen.coeffs.size() >= 1;
    }
    return ok && generators.size() == n * k;
}

size_t Code::maxSize(size_t i) {
    size_t size = 0;
    for(size_t j = i * n; j < (i+1)*n; ++j){
        size = std::max(size, generators[j].coeffs.size());
    }
    return size;
}

bool Code::isSelfOrthogonal(){
    if (k != 1)
        throw new std::logic_error("Not implemented for k > 1");
    Series s;
    for (size_t i = 0; i < n; ++i){
        s = s + generators[i].inverse().conj() * generators[i];
    }
    s.strip();
    if (s.coeffs.size() == 1 && s.coeffs[0] == 0)
        return true;
    return false;
}

bool Code::isOrthogonal(Code& other){
    if (n != other.n)
        throw new std::logic_error("Codes have different n");
    for (size_t i = 0; i < k; ++i){
        for (size_t j = 0; j < other.k; ++j){
            Series s;
            size_t jj = j*n;
            for (size_t ii = i*n; ii < (i+1)*n; ++ii, ++jj){
                s = s + generators[ii].inverse().conj() * other.generators[jj];
            }
            s.strip();
            if (s.coeffs.size() != 1 || s.coeffs[0] != 0)
                return false;
        }
    }
    return true;
}

xlong Code::minDistance(){
    if (!validate()){
        throw std::logic_error("Invalid code");
    }
    std::vector<size_t> lags;
    for (size_t i = 0; i < k; ++i){
        auto l = maxSize(i) - 1;
        if (l < 0)
            throw std::logic_error("Invalid code");
        lags.push_back((size_t)l);
    }

    std::unordered_map<xlong, xlong> d;
    d[0] = INFTY;
    std::priority_queue<std::pair<xlong, xlong>> queue;
    queue.push(std::make_pair<xlong, xlong>(0, 0));
    
    bool firstzero = true;
    while (!queue.empty()){
        xlong v = queue.top().second, curd = -queue.top().first;
        queue.pop();
        
        if (!firstzero && v == 0) { return curd; }
        if (curd > d[v]) continue;
        std::vector<std::vector<gf4>> state(k);
        for (size_t i = 0; i < k; ++i){
            for (size_t j = 0; j < lags[i]; ++j){
                state[i].push_back(gf4((short)(v % GF4_SIZE)));
                v /= GF4_SIZE;
            }
            state[i].push_back(gf4());
        }

        std::vector<gf4> step(k);
        bool allzeros = true;
        for(;;){
            // transition cost and next state
            if (!allzeros || !firstzero){
                xlong nextv(0), t(0);
                xlong pow = 1;
                for (size_t i = 0; i < k; ++i){
                    state[i].back() = step[i];
                    for (size_t j = 1; j < state[i].size(); ++j){
                        nextv = nextv + state[i][j].value * pow;
                        pow *= GF4_SIZE;
                    }
                }
                for (size_t i = 0; i < n; ++i){
                    gf4 s;
                    for (size_t j = 0; j < k; ++j){
                        s = s + generators[j * n + i].dotProduct(state[j]);
                    }
                    if (!(s == 0))
                        t += 1;
                }
                
                auto u = d.find(nextv);
                if (u == d.end()){
                    d[nextv] = curd + t;
                    queue.push(std::make_pair<xlong, xlong>(-(curd + t), xlong(nextv)));
                } else {
                    if (curd + t < u->second){
                        d[nextv] = curd + t;
                        queue.push(std::make_pair<xlong, xlong>(-(curd + t), xlong(nextv)));
                    }
                }
            }

            // next edge
            short carry = 0;
            for (size_t i = 0; i < k; ++i){
                carry = (step[i].value + 1) / GF4_SIZE;
                step[i].value = (step[i].value + 1) % GF4_SIZE;
                if (carry == 0)
                    break;
            }
            if (carry == 1)
                break;
            allzeros = false;
        }
        firstzero = false;
    }

    return -1;
}

std::string Code::toString() {
    std::string s = "";
    for (size_t i = 0; i < generators.size(); ++i){
        s.append(generators[i].toString());
        if (i + 1 != generators.size()){
            if ((i+1) % n == 0)
                s.append("||");
            else 
                s.append("|");
        }
        
    }
    return s;
}

Code Code::findOrthogonal() {
    if (k != 1)
        throw new std::logic_error("Not implemented for k > 1");
    size_t nu = maxSize(0) - 1;
    CodeGenerator code_generator(n, n-k, nu);
    std::shared_ptr<Code> c;
    while ((c = code_generator.next()) != nullptr){
        if (isOrthogonal(*c)){
            return *c;
        }
    }
    return Code(0);
}

CodeGenerator::CodeGenerator(size_t n_, size_t k_, size_t nu_)
: n(n_) 
, k(k_)
, nu(nu_)
, current_degree_split(0)
, current_code(n_ * (k_ + nu_))
{
    degrees = Combinations(nu + k, k).get();
}

std::pair<std::shared_ptr<Code>, bool> CodeGenerator::next_candidate() {
    if (current_degree_split >= degrees->size()){
        current_degree_split = 0;
        bool carry = true;
        size_t i = current_code.size() - 1;
        while (i > 0 && carry){
            carry = (++current_code[i]) == 0;
            --i;
        }
        if (carry)
            return std::make_pair<>(nullptr, true);
    }

    auto code = std::make_shared<Code>(n, k);
    
    auto& v(degrees->at(current_degree_split));
    size_t i = 0;
    for (auto x = v.begin(); x != v.end(); ++x) {
        bool valid = false;
        for (size_t _i = 0; _i < n; ++_i){
            Series s;
            s.coeffs.pop_back();
            for (size_t j = 0; j < *x; ++j, ++i)
                s.coeffs.push_back(current_code[i]);
            valid = valid || s.coeffs.back() != 0;
            code->generators.push_back(s);
        }
        if (!valid) {
            // this is not a valid code, continue the search
            // the degree is less than anticipated
            current_degree_split += 1;
            return std::make_pair<>(nullptr, false);
        }
    }

    current_degree_split += 1;
    return std::make_pair<>(code, true);
}

std::shared_ptr<Code> CodeGenerator::next() {
    std::pair<std::shared_ptr<Code>, bool> c;
    while (!(c = next_candidate()).second);
    return c.first;
}