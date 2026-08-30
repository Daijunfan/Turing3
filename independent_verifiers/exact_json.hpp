#pragma once

#include <cctype>
#include <cstdint>
#include <climits>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fusion_verify {

struct Json {
  enum Kind { Null, Bool, Number, String, Array, Object } kind = Null;
  bool boolean = false;
  long long number = 0;
  std::string string;
  std::vector<Json> array;
  std::unordered_map<std::string, Json> object;

  const Json &at(const std::string &key) const {
    auto it = object.find(key);
    if (it == object.end()) throw std::runtime_error("missing JSON key: " + key);
    return it->second;
  }
};

class JsonParser {
 public:
  explicit JsonParser(std::string text) : text_(std::move(text)) {}
  Json parse() {
    Json value = parseValue();
    whitespace();
    if (pos_ != text_.size()) throw std::runtime_error("trailing JSON data");
    return value;
  }

 private:
  std::string text_;
  size_t pos_ = 0;

  void whitespace() { while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_]))) ++pos_; }
  char take() { if (pos_ >= text_.size()) throw std::runtime_error("unexpected JSON EOF"); return text_[pos_++]; }
  bool consume(char c) { whitespace(); if (pos_ < text_.size() && text_[pos_] == c) { ++pos_; return true; } return false; }

  Json parseValue() {
    whitespace();
    if (pos_ >= text_.size()) throw std::runtime_error("unexpected JSON EOF");
    char c = text_[pos_];
    if (c == '{') return parseObject();
    if (c == '[') return parseArray();
    if (c == '"') { Json j; j.kind = Json::String; j.string = parseString(); return j; }
    if (c == '-' || std::isdigit(static_cast<unsigned char>(c))) return parseNumber();
    if (text_.compare(pos_, 4, "true") == 0) { pos_ += 4; Json j; j.kind = Json::Bool; j.boolean = true; return j; }
    if (text_.compare(pos_, 5, "false") == 0) { pos_ += 5; Json j; j.kind = Json::Bool; return j; }
    if (text_.compare(pos_, 4, "null") == 0) { pos_ += 4; return Json{}; }
    throw std::runtime_error("invalid JSON token at " + std::to_string(pos_));
  }

  Json parseObject() {
    Json j; j.kind = Json::Object; take();
    if (consume('}')) return j;
    do {
      whitespace();
      std::string key = parseString();
      if (!consume(':')) throw std::runtime_error("expected JSON colon");
      j.object.emplace(std::move(key), parseValue());
    } while (consume(','));
    if (!consume('}')) throw std::runtime_error("expected JSON object end");
    return j;
  }

  Json parseArray() {
    Json j; j.kind = Json::Array; take();
    if (consume(']')) return j;
    do { j.array.push_back(parseValue()); } while (consume(','));
    if (!consume(']')) throw std::runtime_error("expected JSON array end");
    return j;
  }

  std::string parseString() {
    if (take() != '"') throw std::runtime_error("expected JSON string");
    std::string out;
    while (true) {
      char c = take();
      if (c == '"') return out;
      if (c != '\\') { out += c; continue; }
      char escaped = take();
      if (escaped == 'u') { for (int i = 0; i < 4; ++i) take(); out += '?'; }
      else if (escaped == 'n') out += '\n';
      else if (escaped == 'r') out += '\r';
      else if (escaped == 't') out += '\t';
      else out += escaped;
    }
  }

  Json parseNumber() {
    size_t start = pos_;
    if (text_[pos_] == '-') ++pos_;
    while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_]))) ++pos_;
    Json j; j.kind = Json::Number; j.number = std::stoll(text_.substr(start, pos_ - start)); return j;
  }
};

inline Json loadJson(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open certificate: " + path);
  return JsonParser(std::string(std::istreambuf_iterator<char>(input), {})).parse();
}

struct Rational {
  long long n = 0, d = 1;
  Rational() = default;
  Rational(long long value) : n(value), d(1) {}
  Rational(long long numerator, long long denominator) : n(numerator), d(denominator) { normalize(); }
  void normalize() {
    if (d == 0) throw std::runtime_error("zero rational denominator");
    if (n == 0) { d = 1; return; }
    if (d < 0) { n = -n; d = -d; }
    long long g = std::gcd(n < 0 ? -n : n, d); n /= g; d /= g;
  }
  bool zero() const { return n == 0; }
};
inline long long narrow(__int128 value) {
  if (value > INT64_MAX || value < INT64_MIN) throw std::runtime_error("exact rational overflow");
  return static_cast<long long>(value);
}
inline Rational operator+(const Rational &a, const Rational &b) {
  long long g=std::gcd(a.d,b.d), ad=a.d/g, bd=b.d/g;
  return Rational(narrow(static_cast<__int128>(a.n)*bd+static_cast<__int128>(b.n)*ad), narrow(static_cast<__int128>(ad)*b.d));
}
inline Rational operator-(const Rational &a, const Rational &b) { return a+Rational(-b.n,b.d); }
inline Rational operator*(const Rational &a, const Rational &b) {
  long long g1=std::gcd(a.n<0?-a.n:a.n,b.d),g2=std::gcd(b.n<0?-b.n:b.n,a.d);
  return Rational(narrow(static_cast<__int128>(a.n/g1)*(b.n/g2)), narrow(static_cast<__int128>(a.d/g2)*(b.d/g1)));
}
inline std::string str(const Rational &a) { return a.d == 1 ? std::to_string(a.n) : std::to_string(a.n)+"/"+std::to_string(a.d); }

struct Coeff { int index; long long n; long long d; };
using Row = std::vector<Coeff>;
struct Product { Row u, v, w; };

inline Row parseRow(const Json &json) {
  Row row;
  for (const Json &triple : json.array) {
    if (triple.array.size() != 3) throw std::runtime_error("sparse coefficient is not [index,numerator,denominator]");
    row.push_back({static_cast<int>(triple.array[0].number), triple.array[1].number, triple.array[2].number});
  }
  return row;
}
inline std::vector<Product> parseProducts(const Json &json) {
  std::vector<Product> products;
  products.reserve(json.array.size());
  for (const Json &item : json.array) products.push_back({parseRow(item.at("U")), parseRow(item.at("V")), parseRow(item.at("W"))});
  return products;
}
inline std::vector<int> intArray(const Json &json) {
  std::vector<int> result; for (const Json &item : json.array) result.push_back(static_cast<int>(item.number)); return result;
}

constexpr long long MOD = 1000003;
inline long long power(long long a, long long e) { long long r=1; for (;e;e>>=1,a=a*a%MOD) if(e&1) r=r*a%MOD; return r; }
inline long long mod(const Coeff &c) { long long n=(c.n%MOD+MOD)%MOD, d=(c.d%MOD+MOD)%MOD; if(!d) throw std::runtime_error("denominator zero modulo screen prime"); return n*power(d,MOD-2)%MOD; }
inline uint64_t key(int a, int b) { return (uint64_t(uint32_t(a))<<32)|uint32_t(b); }
inline std::pair<int,int> unkey(uint64_t value) { return {int(value>>32), int(value&0xffffffffu)}; }

inline size_t verifyStandardMod(const std::vector<Product> &products, int m, int n, int p, bool mutate, std::string &sample) {
  size_t residuals = 0;
  for (int c=0;c<m*p;++c) {
    std::unordered_map<uint64_t,long long> accum;
    for (size_t product=0;product<products.size();++product) {
      for (const Coeff &wc:products[product].w) if(wc.index==c) {
        long long wv=mod(wc);
        for (size_t ui=0;ui<products[product].u.size();++ui) {
          Coeff uc=products[product].u[ui];
          if(mutate && product==0 && ui==0) ++uc.n;
          long long uw=mod(uc)*wv%MOD;
          for(const Coeff &vc:products[product].v) {
            uint64_t k=key(uc.index,vc.index);
            long long value=(accum[k]+uw*mod(vc))%MOD;
            if(value) accum[k]=value; else accum.erase(k);
          }
        }
      }
    }
    int ci=c/p,cj=c%p;
    for(int x=0;x<n;++x) { uint64_t k=key(ci*n+x,x*p+cj); long long value=(accum[k]+MOD-1)%MOD; if(value) accum[k]=value; else accum.erase(k); }
    residuals += accum.size();
    if(sample.empty() && !accum.empty()) { auto [a,b]=unkey(accum.begin()->first); sample="a="+std::to_string(a)+",b="+std::to_string(b)+",c="+std::to_string(c)+",value="+std::to_string(accum.begin()->second); }
  }
  return residuals;
}

inline size_t verifyStandardExact(const std::vector<Product> &products, int m, int n, int p, std::string &sample) {
  size_t residuals=0;
  for(int c=0;c<m*p;++c) {
    std::unordered_map<uint64_t,Rational> accum;
    for(const Product &product:products) for(const Coeff &wc:product.w) if(wc.index==c)
      for(const Coeff &uc:product.u) for(const Coeff &vc:product.v) {
        uint64_t k=key(uc.index,vc.index);
        Rational value=accum[k]+Rational(uc.n,uc.d)*Rational(vc.n,vc.d)*Rational(wc.n,wc.d);
        if(value.zero()) accum.erase(k); else accum[k]=value;
      }
    int ci=c/p,cj=c%p;
    for(int x=0;x<n;++x) { uint64_t k=key(ci*n+x,x*p+cj); Rational value=accum[k]-Rational(1); if(value.zero()) accum.erase(k); else accum[k]=value; }
    residuals += accum.size();
    if(sample.empty() && !accum.empty()) { auto [a,b]=unkey(accum.begin()->first); sample="a="+std::to_string(a)+",b="+std::to_string(b)+",c="+std::to_string(c)+",value="+str(accum.begin()->second); }
  }
  return residuals;
}
}  // namespace fusion_verify
