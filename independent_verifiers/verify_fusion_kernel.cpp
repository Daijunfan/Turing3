#include "exact_json.hpp"
#include <set>

using namespace fusion_verify;

struct Embedding { Row u,v,w; };

size_t verifyKernelMod(const std::vector<Product> &products, const std::vector<Embedding> &slots,
                       const std::vector<int> &ambient, const std::vector<int> &outer,
                       const std::vector<int> &inner, bool mutate, std::string &sample) {
  int on=outer[1],op=outer[2],im=inner[0],in=inner[1],ip=inner[2];
  int pn=ambient[1],pp=ambient[2]; size_t residuals=0;
  for(int c=0;c<ambient[0]*ambient[2];++c) {
    std::unordered_map<uint64_t,long long> accum;
    for(size_t product=0;product<products.size();++product) for(const Coeff &wc:products[product].w) if(wc.index==c)
      for(size_t ui=0;ui<products[product].u.size();++ui) {
        Coeff uc=products[product].u[ui]; if(mutate&&product==0&&ui==0) ++uc.n;
        for(const Coeff &vc:products[product].v) { uint64_t k=key(uc.index,vc.index); long long value=(accum[k]+mod(uc)*mod(vc)%MOD*mod(wc))%MOD; if(value)accum[k]=value;else accum.erase(k); }
      }
    for(const Embedding &slot:slots) for(const Coeff &uo:slot.u) for(const Coeff &vo:slot.v) for(const Coeff &wo:slot.w) {
      int oi=uo.index/on,ok=uo.index%on,ok2=vo.index/op,oj=vo.index%op,oi2=wo.index/op,oj2=wo.index%op;
      long long coeff=mod(uo)*mod(vo)%MOD*mod(wo)%MOD;
      for(int ii=0;ii<im;++ii) for(int ik=0;ik<in;++ik) for(int ij=0;ij<ip;++ij) {
        int cc=(oi2*im+ii)*pp+(oj2*ip+ij); if(cc!=c) continue;
        int a=(oi*im+ii)*pn+(ok*in+ik), b=(ok2*in+ik)*pp+(oj*ip+ij); uint64_t k=key(a,b);
        long long value=(accum[k]+MOD-coeff)%MOD; if(value)accum[k]=value;else accum.erase(k);
      }
    }
    residuals+=accum.size(); if(sample.empty()&&!accum.empty()){auto[a,b]=unkey(accum.begin()->first);sample="a="+std::to_string(a)+",b="+std::to_string(b)+",c="+std::to_string(c)+",value="+std::to_string(accum.begin()->second);}
  }
  return residuals;
}

size_t verifyKernelExact(const std::vector<Product> &products, const std::vector<Embedding> &slots,
                         const std::vector<int> &ambient, const std::vector<int> &outer,
                         const std::vector<int> &inner, std::string &sample) {
  int on=outer[1],op=outer[2],im=inner[0],in=inner[1],ip=inner[2],pn=ambient[1],pp=ambient[2]; size_t residuals=0;
  for(int c=0;c<ambient[0]*ambient[2];++c) {
    std::unordered_map<uint64_t,Rational> accum;
    for(const Product &product:products) for(const Coeff &wc:product.w) if(wc.index==c) for(const Coeff &uc:product.u) for(const Coeff &vc:product.v) {
      uint64_t k=key(uc.index,vc.index); Rational value=accum[k]+Rational(uc.n,uc.d)*Rational(vc.n,vc.d)*Rational(wc.n,wc.d); if(value.zero())accum.erase(k);else accum[k]=value;
    }
    for(const Embedding &slot:slots) for(const Coeff &uo:slot.u) for(const Coeff &vo:slot.v) for(const Coeff &wo:slot.w) {
      int oi=uo.index/on,ok=uo.index%on,ok2=vo.index/op,oj=vo.index%op,oi2=wo.index/op,oj2=wo.index%op; Rational coeff=Rational(uo.n,uo.d)*Rational(vo.n,vo.d)*Rational(wo.n,wo.d);
      for(int ii=0;ii<im;++ii) for(int ik=0;ik<in;++ik) for(int ij=0;ij<ip;++ij) { int cc=(oi2*im+ii)*pp+(oj2*ip+ij);if(cc!=c)continue;int a=(oi*im+ii)*pn+(ok*in+ik),b=(ok2*in+ik)*pp+(oj*ip+ij);uint64_t k=key(a,b);Rational value=accum[k]-coeff;if(value.zero())accum.erase(k);else accum[k]=value; }
    }
    residuals+=accum.size();if(sample.empty()&&!accum.empty()){auto[a,b]=unkey(accum.begin()->first);sample="a="+std::to_string(a)+",b="+std::to_string(b)+",c="+std::to_string(c)+",value="+str(accum.begin()->second);}
  }
  return residuals;
}

int main(int argc,char **argv){
  try{
    if(argc<2||argc>3){std::cerr<<"usage: verify_fusion_kernel CERTIFICATE [--mutate-first]\n";return 2;} bool mutate=argc==3&&std::string(argv[2])=="--mutate-first";
    Json root=loadJson(argv[1]);auto ambient=intArray(root.at("ambient_shape")),outer=intArray(root.at("outer_shape")),inner=intArray(root.at("slot_shape"));auto products=parseProducts(root.at("products"));
    if(products.size()!=static_cast<size_t>(root.at("rank").number))throw std::runtime_error("rank mismatch");std::vector<Embedding> slots;for(const Json &item:root.at("slot_embeddings").array)slots.push_back({parseRow(item.at("U")),parseRow(item.at("V")),parseRow(item.at("W"))});
    std::string sample;size_t modular=verifyKernelMod(products,slots,ambient,outer,inner,mutate,sample);if(modular){std::cout<<"CXX_FUSION_KERNEL_VERIFICATION=FAIL residual_count_mod_1000003="<<modular<<" sample="<<sample<<"\n";return 1;}
    size_t exact=verifyKernelExact(products,slots,ambient,outer,inner,sample);if(exact){std::cout<<"CXX_FUSION_KERNEL_VERIFICATION=FAIL residual_count_exact="<<exact<<" sample="<<sample<<"\n";return 1;}
    long long gain=static_cast<long long>(slots.size())*15-static_cast<long long>(products.size());std::cout<<"CXX_FUSION_KERNEL_VERIFICATION=PASS exact=true rank="<<products.size()<<" fusion_gain="<<gain<<"\n";return 0;
  }catch(const std::exception &e){std::cout<<"CXX_FUSION_KERNEL_VERIFICATION=FAIL error="<<e.what()<<"\n";return 1;}
}
