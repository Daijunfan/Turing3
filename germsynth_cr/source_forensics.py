from __future__ import annotations

import bz2
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def http_headers(url: str) -> dict:
    run = subprocess.run(["curl", "-L", "-sSI", "--connect-timeout", "15", "--max-time", "60", url],
                         text=True, capture_output=True)
    result = {"status": "PASS" if run.returncode == 0 else "FAIL", "error": run.stderr.strip() or None}
    for line in run.stdout.splitlines():
        if ":" not in line:
            if line.startswith("HTTP/"):
                result["http_status"] = line.strip()
            continue
        key, value = line.split(":", 1)
        key = key.lower()
        if key in ("etag", "last-modified", "content-length", "content-type"):
            result[key.replace("-", "_")] = value.strip()
    return result


def file_record(url: str, path: str | Path, license_status: str) -> dict:
    path = Path(path)
    try:
        display_path = str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        display_path = str(path)
    compressed = path.read_bytes()
    try:
        expanded = bz2.decompress(compressed)
        bzip2 = "PASS"
    except OSError as error:
        expanded, bzip2 = b"", f"FAIL: {error}"
    return {"url": url, "path": display_path, "downloaded_at": datetime.now(timezone.utc).isoformat(),
            **http_headers(url), "compressed_sha256": sha256_bytes(compressed),
            "decompressed_sha256": sha256_bytes(expanded) if expanded else None,
            "compressed_bytes": len(compressed), "decompressed_bytes": len(expanded),
            "bzip2_integrity": bzip2, "license": license_status}


def plain_file_record(url: str, path: str | Path, license_status: str) -> dict:
    path = Path(path); data = path.read_bytes()
    try:
        display_path = str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        display_path = str(path)
    return {"url": url, "path": display_path, "downloaded_at": datetime.now(timezone.utc).isoformat(),
            **http_headers(url), "sha256": sha256_bytes(data), "bytes": len(data),
            "license": license_status}


def history_probe(url: str) -> dict:
    endpoint = ("https://web.archive.org/cdx/search/cdx?url=" + url +
                "&output=json&fl=timestamp,original,statuscode,digest,length&filter=statuscode:200&collapse=digest")
    run = subprocess.run(["curl", "-L", "-sS", "--connect-timeout", "10", "--max-time", "20", endpoint],
                         text=True, capture_output=True)
    if run.returncode:
        return {"status": "FAIL", "service": "Internet Archive CDX", "error": run.stderr.strip()}
    try:
        return {"status": "PASS", "service": "Internet Archive CDX", "captures": json.loads(run.stdout)}
    except json.JSONDecodeError as error:
        return {"status": "FAIL", "service": "Internet Archive CDX", "error": str(error),
                "response_prefix": run.stdout[:500]}


def write_lock(records: list[dict], repositories: list[dict], history: list[dict], path="sources.lock.json"):
    lock = {"format": "germsynth-cr-sources-lock-v1", "generated_at": datetime.now(timezone.utc).isoformat(),
            "files": records, "repositories": repositories, "history_queries": history}
    Path(path).write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    return lock
