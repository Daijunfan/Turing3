#include "exact_json.hpp"
using namespace fusion_verify;
int main(int argc,char**argv){
 try{
  if(argc!=2){std::cerr<<"usage: verify_residual_completion CERTIFICATE\n";return 2;}
  Json root=loadJson(argv[1]);
  if(root.at("status").string!="PASS"){
   std::cout<<"CXX_RESIDUAL_COMPLETION=FAIL certificate_status="<<root.at("status").string
            <<" lower_bound="<<root.at("lower_bound").number<<"\n";return 1;
  }
  if(root.at("factors").array.size()!=static_cast<size_t>(root.at("rank").number))throw std::runtime_error("completion rank mismatch");
  std::cout<<"CXX_RESIDUAL_COMPLETION=FAIL error=PASS completion schema requires embedded residual coordinates\n";return 1;
 }catch(const std::exception&e){std::cout<<"CXX_RESIDUAL_COMPLETION=FAIL error="<<e.what()<<"\n";return 1;}
}
