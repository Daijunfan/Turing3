#include "exact_json.hpp"
using namespace fusion_verify;
int main(int argc,char**argv){
 try{
  if(argc!=2){std::cerr<<"usage: verify_novel_candidate MANIFEST_OR_CERTIFICATE\n";return 2;}
  Json root=loadJson(argv[1]);
  auto it=root.object.find("candidates");
  if(it!=root.object.end()&&it->second.array.empty()){std::cout<<"CXX_NOVEL_CANDIDATE=FAIL no_candidate\n";return 1;}
  std::cout<<"CXX_NOVEL_CANDIDATE=FAIL unsupported_or_unverified_candidate\n";return 1;
 }catch(const std::exception&e){std::cout<<"CXX_NOVEL_CANDIDATE=FAIL error="<<e.what()<<"\n";return 1;}
}
