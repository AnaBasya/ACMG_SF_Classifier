# cooccurrence_helper.py
# Place this file next to acmg_sf_classifier.py and import functions where needed.

import time
import json
import requests
import os
from typing import List, Tuple, Dict, Optional

GRAPHQL_ENDPOINT_DEFAULT = "https://gnomad.broadinstitute.org/api"
DEFAULT_DATASET = "gnomad_r2_1"
SESSION = requests.Session()
SESSION.headers.update({"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "cooccurrence-fetcher/1.0"})

# A compact GraphQL template that works for common gnomAD deployments.
# This template may need adjustment if gnomAD schema differs; we attempt fallback handling.
GRAPHQL_QUERY_VARIANT_CO = """
query PairCooccurrence($variantId1: String!, $variantId2: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId1, dataset: $dataset) {
    variantId
    cooccurrence(targetVariantId: $variantId2) {
      n_carriers
      n_cis
      n_trans
      observations_count
      # Note: the web API may not expose raw sample-level list; we rely on aggregated counts.
    }
  }
}
"""

GRAPHQL_QUERY_CO_FALLBACK = """
query PairCooccurrenceFallback($variantId1: String!, $variantId2: String!, $dataset: DatasetId!) {
  cooccurrence(variant1: $variantId1, variant2: $variantId2, dataset: $dataset) {
    n_carriers
    n_cis
    n_trans
    evidence
  }
}
"""

def _graphql_query(endpoint: str, query: str, variables: dict, timeout: int = 30) -> dict:
    payload = {"query": query, "variables": variables}
    for attempt in range(1, 5):
        try:
            r = SESSION.post(endpoint, json=payload, timeout=timeout)
            if r.status_code == 429:
                wait = int(r.headers.get("Retry-After", 5))
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except requests.RequestException as ex:
            # exponential backoff
            sleep_for = min(60, 2 ** attempt)
            time.sleep(sleep_for)
    raise RuntimeError("GraphQL query failed after retries")

def query_pair_gnomad(vid1: str, vid2: str, dataset: str = DEFAULT_DATASET,
                      endpoint: str = GRAPHQL_ENDPOINT_DEFAULT, use_fallback: bool = False) -> Optional[Dict]:
    """
    Query gnomAD for cooccurrence of a single pair. Returns a dict or None.
    vid format: 'chr1:12345:A:T'
    """
    vars = {"variantId1": vid1, "variantId2": vid2, "dataset": dataset}
    # try primary template first
    try:
        q = GRAPHQL_QUERY_CO_FALLBACK if use_fallback else GRAPHQL_QUERY_VARIANT_CO
        resp = _graphql_query(endpoint, q, vars)
        data = resp.get("data")
        if not data:
            return None
        # try to extract cooccurrence info from common places:
        if data.get("variant") and isinstance(data["variant"], dict):
            co = data["variant"].get("cooccurrence")
            if co is not None:
                return co
            # fallback - sometimes variant returns other shapes
            return data["variant"]
        if data.get("cooccurrence"):
            return data.get("cooccurrence")
        return data
    except Exception:
        # second attempt with fallback template
        if not use_fallback:
            return query_pair_gnomad(vid1, vid2, dataset=dataset, endpoint=endpoint, use_fallback=True)
        return None

def fetch_cooccurrence_for_pairs(pairs: List[Tuple[str,str]],
                                 dataset: str = DEFAULT_DATASET,
                                 endpoint: str = GRAPHQL_ENDPOINT_DEFAULT,
                                 delay: float = 0.2,
                                 max_queries: Optional[int] = None) -> Dict[Tuple[str,str], Dict]:
    """
    Fetch cooccurrence info for the given list of unordered pairs.
    Returns dict mapping (v1,v2) -> cooccurrence_info(dict) or {'error':...}.
    Caches results during this run (local only).
    """
    results: Dict[Tuple[str,str], Dict] = {}
    done = 0
    for v1, v2 in pairs:
        key = (v1, v2)
        if max_queries and done >= max_queries:
            break
        try:
            info = query_pair_gnomad(v1, v2, dataset=dataset, endpoint=endpoint)
            if info is None:
                results[key] = {"n_cis": None, "n_trans": None, "n_carriers": None, "raw": None}
            else:
                # normalize possible fields
                n_cis = info.get("n_cis") if isinstance(info, dict) else None
                n_trans = info.get("n_trans") if isinstance(info, dict) else None
                n_carriers = info.get("n_carriers") if isinstance(info, dict) else None
                results[key] = {"n_cis": n_cis, "n_trans": n_trans, "n_carriers": n_carriers, "raw": info}
        except Exception as e:
            results[key] = {"error": str(e)}
        done += 1
        time.sleep(delay)
    return results
