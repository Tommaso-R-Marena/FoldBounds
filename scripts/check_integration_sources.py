"""Fail-closed, dependency-free provenance/structure sanity for the integrated research tree."""
from __future__ import annotations
import json
import re
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    data = json.loads((root / "research/integration_sources.json").read_text(encoding="utf-8"))
    tips = data["selected_tips"]
    assert data["schema_version"] == 1
    assert len(tips) == 20, "Expected twenty selected tip branches"
    assert len({item["tip_pr"] for item in tips}) == len(tips)
    assert data["covered_prs"] == sorted(set(sum((item["dependency_chain"] for item in tips), [])))
    assert len(data["covered_prs"]) == 61, "Source PR inventory changed; refresh manifest deliberately"
    assert set(data["excluded_obsolete_prs"]) == {9}
    for item in tips:
        assert re.fullmatch(r"[0-9a-f]{40}", item["head_sha"]), item
        assert item["dependency_chain"][-1] == item["tip_pr"], item
        path = root / item["anchor"]
        assert path.is_file(), f"Missing integrated source for PR #{item['tip_pr']}: {path}"
        assert path.resolve().is_relative_to(root.resolve()), f"Unsafe anchor path: {path}"
    print(f"PASS: {len(tips)} selected tips, {len(data['covered_prs'])} PR dependencies, source anchors present")


if __name__ == "__main__":
    main()
