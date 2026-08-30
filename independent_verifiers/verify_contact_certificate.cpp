#include "exact_json.hpp"
using namespace fusion_verify;
int main(int argc,char**argv){
 try{
  if(argc!=2){std::cerr<<"usage: verify_contact_certificate PARENT_STYLE_CERTIFICATE\n";return 2;}
  Json root=loadJson(argv[1]);auto shape=intArray(root.at("shape"));auto products=parseProducts(root.at("products"));
  if(products.size()!=static_cast<size_t>(root.at("rank").number))throw std::runtime_error("rank mismatch");
  std::string sample;size_t mod=verifyStandardMod(products,shape[0],shape[1],shape[2],false,sample);
  if(mod){std::cout<<"CXX_CONTACT_CERTIFICATE=FAIL residual_mod="<<mod<<" sample="<<sample<<"\n";return 1;}
  size_t exact=verifyStandardExact(products,shape[0],shape[1],shape[2],sample);
  if(exact){std::cout<<"CXX_CONTACT_CERTIFICATE=FAIL residual_exact="<<exact<<" sample="<<sample<<"\n";return 1;}
  std::cout<<"CXX_CONTACT_CERTIFICATE=PASS exact=true rank="<<products.size()<<"\n";return 0;
 }catch(const std::exception&e){std::cout<<"CXX_CONTACT_CERTIFICATE=FAIL error="<<e.what()<<"\n";return 1;}
}
