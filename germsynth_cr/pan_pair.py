from __future__ import annotations

from fractions import Fraction


def _put(row, index, value):
    row[index] = row.get(index, Fraction()) + Fraction(value)
    if not row[index]:
        row.pop(index)


def build(n: int, r: int, p: int):
    du, dv, dw = n*r+p*n, r*p+n*r, n*p+p*r
    terms = []
    a1=lambda i,j:i*r+j
    a2=lambda k,i:n*r+k*n+i
    b1=lambda j,k:j*p+k
    b2=lambda i,j:r*p+i*r+j
    c1=lambda i,k:i*p+k
    c2=lambda k,j:n*p+k*r+j
    for i in range(n):
        for j in range(r):
            for k in range(p):
                U,V,W={},{},{}
                _put(U,a1(i,j),1);_put(U,a2(k,i),1)
                _put(V,b1(j,k),1);_put(V,b2(i,j),1)
                _put(W,c1(i,k),1);_put(W,c2(k,j),1)
                terms.append((U,V,W))
    for i in range(n):
        for k in range(p):
            U,V,W={a2(k,i):Fraction(-1)},{},{c1(i,k):Fraction(1)}
            for j in range(r):_put(V,b1(j,k),1);_put(V,b2(i,j),1)
            terms.append((U,V,W))
    for i in range(n):
        for j in range(r):
            U,V,W={a1(i,j):Fraction(-1)},{b2(i,j):Fraction(1)},{}
            for k in range(p):_put(W,c1(i,k),1);_put(W,c2(k,j),1)
            terms.append((U,V,W))
    for j in range(r):
        for k in range(p):
            U,V,W={}, {b1(j,k):Fraction(1)}, {c2(k,j):Fraction(1)}
            for i in range(n):_put(U,a1(i,j),-1);_put(U,a2(k,i),-1)
            terms.append((U,V,W))
    return {"dimensions": [du,dv,dw], "terms": terms, "shape_pair": [[n,r,p],[p,n,r]],
            "rank": len(terms), "separate_naive_rank": 2*n*r*p,
            "gain_vs_naive": 2*n*r*p-len(terms), "mixed_product_count": n*r*p}


def verify(instance) -> dict:
    actual = {}
    for U,V,W in instance["terms"]:
        for a,uc in U.items():
            for b,vc in V.items():
                for c,wc in W.items():
                    key=(a,b,c);value=actual.get(key,Fraction())+uc*vc*wc
                    if value:actual[key]=value
                    else:actual.pop(key,None)
    n,r,p=instance["shape_pair"][0]
    expected={}
    for i in range(n):
        for j in range(r):
            for k in range(p):expected[(i*r+j,j*p+k,i*p+k)]=Fraction(1)
    for k in range(p):
        for i in range(n):
            for j in range(r):expected[(n*r+k*n+i,r*p+i*r+j,n*p+k*r+j)]=Fraction(1)
    residual={key:actual.get(key,Fraction())-value for key,value in expected.items()}
    for key,value in actual.items():
        if key not in expected:residual[key]=value
    residual={key:value for key,value in residual.items() if value}
    return {"status":"PASS" if not residual else "FAIL","rank":instance["rank"],
            "residual_count":len(residual),"residual_sample":[(*key,str(value)) for key,value in list(residual.items())[:16]],
            "mixed_product_count":instance["mixed_product_count"],"gain_vs_naive":instance["gain_vs_naive"]}
