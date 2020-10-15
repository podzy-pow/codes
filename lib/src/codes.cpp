#include "codes.h"

#define LOG(msg) \
    std::cerr << __FILE__ << "(" << __LINE__ << "): " << msg << std::endl 

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

size_t Code::maxSize(){
    size_t size = 0;
    for(size_t j = 0; j < k*n; ++j){
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

Code Code::findOrthogonalOld() {
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

Code Code::findOrthogonal(){
    // set degrees of each generator in the orthogonal system
    size_t other_k = n - k;
    size_t other_degree = maxSize(0) - 1;
    for (size_t j = 1; j < k; ++j){
        other_degree += maxSize(j) - 1;
    }
    std::vector<size_t> other_degrees;
    size_t avg = other_degree / other_k;
    for (size_t j = 0; j < other_k; ++j){
        other_degrees.push_back(avg);
        other_degree -= avg;
    }
    for (size_t j = 0; j < other_degree; ++j){
        other_degrees[j] += 1;
    }

    for (size_t j = 0; j < other_k; ++j)
        std::cerr<< other_degrees[j] << " ";
    std::cerr << std::endl;


    // define index of each coefficient of the orthogonal system 
    // in system of linear equations as shift + index in generator
    std::vector<size_t> coeffs_shifts;
    std::cerr << "other_k: " << other_k << " " << "n: " << n << std::endl;
    size_t M = 0;
    for (size_t i = 0; i < other_k; ++i){
        for (size_t j = 0; j < n; ++j){
            coeffs_shifts.push_back(M);
            M += other_degrees[i] + 1;
        }
    }
    

    // Preparing system of linear equations
    std::vector<std::vector<gf4>> linear_system;

    size_t N = 0;
    for (size_t i = 0; i < k; ++i){
        size_t size = maxSize(i) - 1;
        for (size_t j = 0; j < other_k; ++j){
            for (size_t z = 0; z <= size + other_degrees[j]; ++z){
                std::vector<gf4> v;
                for (size_t zz = 0; zz < M; ++zz){
                    v.push_back(gf4());
                }
                linear_system.push_back(v);
            }
            for (size_t m = 0; m < n; ++m){
                Series& this_generator = generators[i * n + m];
                size_t other_generator_index = j * n + m;
                for (size_t ii = 0; ii <= size; ++ii){
                    size_t p = size - ii;
                    for (size_t jj = 0; jj <= other_degrees[j]; ++jj){
                        linear_system[N + p + jj][coeffs_shifts[other_generator_index] + jj] = this_generator.at(ii).conj();
                    }
                }
            }
            N += size + other_degrees[j] + 1;
        }
    }

    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < M; ++j){
            std::cerr << linear_system[i][j].toString() << " ";
        }
        std::cerr << std::endl;
    }

    // Upper-triangle
    size_t most_left = 0;
    for (size_t i = 0; i < N; ++i){
        size_t j;
        for(;;){
            j = i;
            while (j < N){
                if (linear_system[j][most_left] != 0)
                    break;
                ++j;
            }
            if (j < N)
                break;
            ++most_left;
            if (most_left == M)
                break;
        }
        if (most_left == M)
            break;
        if (i != j)
            std::swap(linear_system[i], linear_system[j]);
        if (linear_system[i][most_left] != 1){
            gf4 c = linear_system[i][most_left].conj();
            for (size_t j = most_left; j < M; ++j)
                linear_system[i][j] = linear_system[i][j] * c;
        }
        
        for (size_t j = i + 1; j < N; ++j){
            if (linear_system[j][most_left] != 0){
                gf4 c = linear_system[j][most_left];
                for (size_t jj = most_left; jj < M; ++jj)
                    linear_system[j][jj] = linear_system[j][jj] + linear_system[i][jj] * c;
            }
        }

        ++most_left;
    }

    std::cerr << "After getting to upper-triangle" << std::endl;
    for (size_t i = 0; i < N; ++i){
        for (size_t j = 0; j < M; ++j){
            std::cerr << linear_system[i][j].toString() << " ";
        }
        std::cerr << std::endl;
    }

    // Assign values. What is the right way?
    std::vector<gf4> values;
    std::vector<bool> assigned;
    char Z = 1;
    for (size_t i = 0; i < M; ++i){
        values.push_back(gf4());
        assigned.push_back(false);
    }
    for (size_t ii = 0; ii < N; ++ii){
        size_t i = N - ii - 1;
        size_t most_left = 0;
        while (most_left < M && linear_system[i][most_left] == 0){
            ++most_left;
        }
        if (most_left == M){
            continue;
        }
        gf4 v = 0;
        for (size_t j = most_left + 1; j < M; ++j){
            if (!assigned[j]){
                // Here we need to make sure 
                // that generators are linearly independent(!)
                // This is not a good way.
                values[j] = gf4(Z);
                Z = Z + 1;
                if (Z == 4)
                    Z = 0;
                assigned[j] = true;
            }
            v = v + linear_system[i][j] * values[j];
        }
        values[most_left] = v;
        assigned[most_left] = true;
    }

    std::cerr << "Found vector" << std::endl;
    for (size_t i = 0; i < M; ++i)
        std::cerr << values[i].toString() << " ";
    std::cerr << std::endl;

    std::cerr << "Verification" << std::endl;
    for (size_t i = 0; i < N; ++i){
        gf4 v = 0;
        for (size_t j = 0; j < M; ++j)
            v = v + values[j] * linear_system[i][j];
        std::cerr << "i: " << i << ", v: " << v.toString() << std::endl;
    }
    
    std::vector<Series> series;
    for (size_t i = 0; i < other_k; ++i){
        for (size_t j = 0; j < n; ++j){
            std::vector<gf4> gen_coeffs;
            size_t jj = coeffs_shifts[i * n + j];
            for (size_t ii = 0; ii <= other_degrees[i]; ++ii, ++jj){
                gen_coeffs.push_back(values[jj]);
            }
            series.push_back(Series(gen_coeffs));
        }
    }
    
    return Code(series, n, other_k);
}


