from pathlib import Path
import csv
import hashlib


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
EXPECTED = {
    "tables/proteins.csv": 1996,
    "tables/loop_annotations.csv": 3383,
    "tables/loop_flexibility_public_summary.csv": 3383,
    "tables/protein_flexibility_summary.csv": 1996,
}
FORBIDDEN_EXTENSIONS = {
    ".doc",
    ".docx",
    ".jpeg",
    ".jpg",
    ".pdf",
    ".png",
    ".ppt",
    ".pptx",
    ".svg",
    ".tif",
    ".tiff",
    ".xls",
    ".xlsx",
}


def count_rows(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        next(reader)
        return sum(1 for _ in reader)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


failed = False

for rel, expected in EXPECTED.items():
    path = ROOT / rel
    observed = count_rows(path)
    status = "OK" if observed == expected else "FAIL"
    print(f"{status}\ttable\t{rel}\tobserved={observed}\texpected={expected}")
    failed = failed or observed != expected

manifest_path = REPO_ROOT / "MANIFEST_SHA256.csv"
with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
    manifest_rows = list(csv.DictReader(handle))

listed = {row["path"]: row for row in manifest_rows}
actual = {
    path.relative_to(REPO_ROOT).as_posix(): path
    for path in REPO_ROOT.rglob("*")
    if path.is_file()
    and path != manifest_path
    and ".git" not in path.relative_to(REPO_ROOT).parts
}

missing = sorted(set(listed) - set(actual))
unlisted = sorted(set(actual) - set(listed))
bad_hashes = []
bad_sizes = []
for rel, row in listed.items():
    path = actual.get(rel)
    if path is None:
        continue
    if str(path.stat().st_size) != str(row["bytes"]):
        bad_sizes.append(rel)
    if sha256(path) != row["sha256"].lower():
        bad_hashes.append(rel)

manifest_ok = not missing and not unlisted and not bad_hashes and not bad_sizes
print(
    f"{'OK' if manifest_ok else 'FAIL'}\tmanifest\tentries={len(manifest_rows)} "
    f"missing={len(missing)} unlisted={len(unlisted)} "
    f"bad_hashes={len(bad_hashes)} bad_sizes={len(bad_sizes)}"
)
for label, values in (
    ("missing", missing),
    ("unlisted", unlisted),
    ("bad_hash", bad_hashes),
    ("bad_size", bad_sizes),
):
    for value in values:
        print(f"FAIL\tmanifest_{label}\t{value}")
failed = failed or not manifest_ok

for nested_manifest in sorted(
    (path for path in REPO_ROOT.rglob("MANIFEST_SHA256.csv") if path != manifest_path),
    key=lambda path: path.as_posix(),
):
    with nested_manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        nested_rows = list(csv.DictReader(handle))
    nested_listed = {row["path"]: row for row in nested_rows}
    nested_actual = {
        path.relative_to(nested_manifest.parent).as_posix(): path
        for path in nested_manifest.parent.rglob("*")
        if path.is_file() and path != nested_manifest
    }
    nested_missing = sorted(set(nested_listed) - set(nested_actual))
    nested_unlisted = sorted(set(nested_actual) - set(nested_listed))
    nested_bad = []
    for rel, row in nested_listed.items():
        path = nested_actual.get(rel)
        if path is None:
            continue
        if (
            str(path.stat().st_size) != str(row["bytes"])
            or sha256(path) != row["sha256"].lower()
        ):
            nested_bad.append(rel)
    nested_ok = not nested_missing and not nested_unlisted and not nested_bad
    rel_manifest = nested_manifest.relative_to(REPO_ROOT).as_posix()
    print(
        f"{'OK' if nested_ok else 'FAIL'}\tmanifest\t{rel_manifest} "
        f"entries={len(nested_rows)} missing={len(nested_missing)} "
        f"unlisted={len(nested_unlisted)} bad={len(nested_bad)}"
    )
    failed = failed or not nested_ok

forbidden = sorted(
    rel for rel, path in actual.items() if path.suffix.lower() in FORBIDDEN_EXTENSIONS
)
for rel in forbidden:
    print(f"FAIL\tforbidden_file_type\t{rel}")
failed = failed or bool(forbidden)

raise SystemExit(1 if failed else 0)
