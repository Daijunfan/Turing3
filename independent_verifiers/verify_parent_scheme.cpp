#include "exact_json.hpp"

using namespace fusion_verify;

int main(int argc, char **argv) {
  try {
    if (argc < 2 || argc > 3) { std::cerr << "usage: verify_parent_scheme CERTIFICATE [--mutate-first]\n"; return 2; }
    bool mutate = argc == 3 && std::string(argv[2]) == "--mutate-first";
    Json root = loadJson(argv[1]);
    std::vector<int> shape = intArray(root.at("shape"));
    std::vector<Product> products = parseProducts(root.at("products"));
    if (shape.size()!=3 || products.size()!=static_cast<size_t>(root.at("rank").number)) throw std::runtime_error("rank or shape mismatch");
    std::string sample;
    size_t modular = verifyStandardMod(products,shape[0],shape[1],shape[2],mutate,sample);
    if(modular) { std::cout << "CXX_PARENT_INDEPENDENT_VERIFICATION=FAIL residual_count_mod_1000003=" << modular << " sample=" << sample << "\n"; return 1; }
    size_t exact = verifyStandardExact(products,shape[0],shape[1],shape[2],sample);
    if(exact) { std::cout << "CXX_PARENT_INDEPENDENT_VERIFICATION=FAIL residual_count_exact=" << exact << " sample=" << sample << "\n"; return 1; }
    std::cout << "CXX_PARENT_INDEPENDENT_VERIFICATION=PASS exact=true rank=" << products.size() << "\n";
    return 0;
  } catch(const std::exception &e) { std::cout << "CXX_PARENT_INDEPENDENT_VERIFICATION=FAIL error=" << e.what() << "\n"; return 1; }
}
