#include "exact_json.hpp"
using namespace fusion_verify;
int main(int argc,char**argv){
 try{
  if(argc!=2){std::cerr<<"usage: verify_8x27x30_3744 CERTIFICATE\n";return 2;}
  Json root=loadJson(argv[1]);
  if(root.at("status").string!="PASS"){std::cout<<"CXX_RANK_3744_VERIFICATION=FAIL certificate_status="<<root.at("status").string<<"\n";return 1;}
  auto shape=intArray(root.at("shape"));auto products=parseProducts(root.at("products"));
  if(shape!=std::vector<int>({8,27,30})||products.size()!=3744)throw std::runtime_error("expected <8,27,30>:3744");
  std::string sample;size_t mod=verifyStandardMod(products,8,27,30,false,sample);
  if(mod){std::cout<<"CXX_RANK_3744_VERIFICATION=FAIL residual_mod="<<mod<<" sample="<<sample<<"\n";return 1;}
  size_t exact=verifyStandardExact(products,8,27,30,sample);
  if(exact){std::cout<<"CXX_RANK_3744_VERIFICATION=FAIL residual_exact="<<exact<<" sample="<<sample<<"\n";return 1;}
  std::cout<<"CXX_RANK_3744_VERIFICATION=PASS exact=true rank=3744\n";return 0;
 }catch(const std::exception&e){std::cout<<"CXX_RANK_3744_VERIFICATION=FAIL error="<<e.what()<<"\n";return 1;}
}
