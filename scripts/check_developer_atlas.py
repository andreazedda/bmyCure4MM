"""Check the developer source inventory without importing apps or opening data."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    inventory = json.loads((root / "docs/development/SOURCE_INVENTORY.json").read_text())
    source_sha = inventory["source_sha"]
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha):
        raise ValueError("Inventory requires an exact baseline source SHA")
    count = 0
    for module in inventory["modules"]:
        relative = Path(module["path"])
        source = (root / relative).resolve()
        if not source.is_relative_to(root) or source.suffix != ".py" or not source.is_file():
            raise ValueError(f"Missing or invalid source path: {relative}")
        tree = ast.parse(source.read_text())
        symbols = {
            node.name: node
            for node in tree.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
        }
        for entry in module["symbols"]:
            node = symbols.get(entry["name"])
            if node is None:
                raise ValueError(f"Missing inventoried symbol: {relative}:{entry['name']}")
            signature = node.name
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                signature += "(" + ast.unparse(node.args) + ")"
            if signature != entry["signature"]:
                raise ValueError(f"Signature changed; review atlas: {relative}:{node.name}")
            expected = (
                f"https://github.com/andreazedda/bmyCure4MM/blob/{source_sha}/"
                f"{relative.as_posix()}#L{entry['line']}"
            )
            if entry["source_url"] != expected or entry["line"] < 1:
                raise ValueError(f"Invalid baseline source anchor: {relative}:{node.name}")
            if isinstance(node, ast.ClassDef):
                members = {
                    item.name: item
                    for item in node.body
                    if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef))
                }
                for method in entry.get("methods", []):
                    member = members.get(method["name"])
                    if member is None or (
                        member.name + "(" + ast.unparse(member.args) + ")" != method["signature"]
                    ):
                        raise ValueError(
                            f"Method changed; review atlas: {relative}:{node.name}.{method['name']}"
                        )
                fields = {
                    item.targets[0].id: item.value
                    for item in node.body
                    if isinstance(item, ast.Assign)
                    and len(item.targets) == 1
                    and isinstance(item.targets[0], ast.Name)
                    and isinstance(item.value, ast.Call)
                }
                for field in entry.get("declared_fields", []):
                    declaration = fields.get(field["name"])
                    expected_type = ast.parse(field["declaration"], mode="eval").body
                    if declaration is None or ast.unparse(declaration.func) != ast.unparse(
                        expected_type.func
                    ):
                        raise ValueError(
                            f"Field type changed; review atlas: {relative}:{node.name}.{field['name']}"
                        )
            count += 1
    print(f"DEVELOPER_ATLAS=PASS modules={len(inventory['modules'])} symbols={count}")
    print("Static paths/signatures/field types only; baseline anchors are versioned, not live line numbers.")


if __name__ == "__main__":
    main()
