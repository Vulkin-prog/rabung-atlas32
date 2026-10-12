#!/usr/bin/env python3
"""Check actual series editions and the planned chronological citation graph.

The registry records author-supplied public editions, not automatic publication.
It makes no network call and does not infer public availability from a commit.
"""
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def normalized(text):
    text = re.sub(r"\\[A-Za-z]+", "", text).lower()
    return re.sub(r"[^a-z0-9]", "", text)


def main():
    registry = json.loads((ROOT / "docs/PUBLICATION_SEQUENCE.json").read_text())
    manuscripts = {m["id"]: m for m in registry["manuscripts"]}
    editions = {p["id"]: p for p in registry["publications"]}
    require(len(manuscripts) == 4, "Expected the four-article portfolio")
    require(len(editions) == len(registry["publications"]), "Repeated edition ID")
    require(sorted(m["release_order"] for m in manuscripts.values()) == [1, 2, 3, 4],
            "Ambiguous release order")
    for insertion in registry["planned_insertions_after_publication"]:
        require(manuscripts[insertion["after"]]["release_order"] <
                manuscripts[insertion["in"]]["release_order"],
                "Planned forward citation")
    checked = 0
    for m in manuscripts.values():
        source = (ROOT / m["source"]).read_text()
        require(m.get("edition_version", m.get("draft_version")) in source, "Edition version differs: " + m["id"])
        # TeX comments can retain editorial instructions; they are not printed references.
        printed = re.sub(r"(?<!\\)%[^\n]*", "", source)
        marker = r"\\bibitem(?:\[[\s\S]*?\])?\{([^}]+)\}"
        matches = list(re.finditer(marker, printed))
        items = {}
        for i, match in enumerate(matches):
            end = matches[i+1].start() if i+1 < len(matches) else len(printed)
            items[match.group(1)] = printed[match.end():end]
        cited = {key.strip() for match in re.finditer(
            r"\\cite[A-Za-z]*\s*(?:\[[^\]]*\]\s*)*\{([^}]+)\}", printed)
                 for key in match.group(1).split(",")}
        require(cited <= items.keys(), "Undefined citation keys: " + m["id"])
        known = {r["key"]: r["publication"] for r in m["series_references"]}
        require(len(known) == len(m["series_references"]), "Repeated series citation key")
        found = set()
        for key, body in items.items():
            related = [p for p in manuscripts.values()
                       if normalized(p["title"]) in normalized(body)]
            if not related:
                require(key not in known, "Registered series citation has wrong title")
                continue
            require(len(related) == 1 and key in known,
                    "Unregistered series bibliography entry: " + m["id"] + ":" + key)
            edition = editions[known[key]]
            require(edition["paper"] == related[0]["id"] and edition["status"] == "published",
                    "Series entry is not a published edition")
            require(edition["release_order"] < m["release_order"], "Forward series reference")
            require(edition["version"] in body and edition["url"] in body,
                    "Exact edition/version locator missing")
            require(key in cited, "Unused series bibliography entry")
            found.add(key)
            checked += 1
        require(found == known.keys(), "Series reference registry differs: " + m["id"])
    print(json.dumps({"status": "PASS", "manuscripts": len(manuscripts),
                      "existing_public_editions": len(editions),
                      "series_references_checked": checked, "forward_references": 0,
                      "scope": "declared public editions and bibliography; no automatic publication"}))


if __name__ == "__main__":
    main()
