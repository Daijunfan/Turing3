from __future__ import annotations

import hashlib
import json


def canonical_sparse_products(products) -> str:
    encoded = []
    for U, V, W in products:
        factors = [tuple((index, str(value)) for index, value in sorted(row.items())) for row in (U, V, W)]
        encoded.append(tuple(factors))
    payload = json.dumps(sorted(encoded), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()
