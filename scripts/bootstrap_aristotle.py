#!/usr/bin/env python3
"""One-time, allowlisted import from a checksum-pinned Aristotle archive.

No extracted file is executed or imported as Python. The expected archive
format is .tar.xz containing precisely the pinned Lean package sources.
"""
from pathlib import Path
import base64
import hashlib
import io
import lzma
import tarfile

ROOT=Path(".")
PARTS=[ROOT/".aristotle-staging"/("part"+str(i)+".b64") for i in range(5)]
ARCHIVE_SHA256="eae16a0982885a213b8a003bb6bb534ccfb40d57eef6e9f189e711211cf57735"
SOURCES={
    "lean-toolchain", "lakefile.toml", "lake-manifest.json", "FoldBounds.lean",
    "FoldBounds/Basic.lean", "FoldBounds/PartitionBounds.lean",
    "FoldBounds/Probability.lean", "FoldBounds/FreeEnergy.lean",
    "FoldBounds/Compositional.lean", "FoldBounds/Symmetry.lean",
    "FoldBounds/Certificate.lean", "FoldBounds/AxiomAudit.lean",
    "FoldBounds/Optimality.lean", "FoldBounds/Variational.lean",
}
def main():
    payload="".join(p.read_text(encoding="ascii") for p in PARTS)
    archive=base64.b64decode(payload,validate=True)
    actual=hashlib.sha256(archive).hexdigest()
    if actual != ARCHIVE_SHA256:
        raise RuntimeError("Aristotle source transfer checksum failed: "+actual)
    data=lzma.decompress(archive)
    with tarfile.open(fileobj=io.BytesIO(data),mode="r:") as tar:
        members=[m for m in tar.getmembers() if m.isfile()]
        names={m.name for m in members}
        if names!=SOURCES or any(not m.isfile() for m in members):
            raise RuntimeError("Unexpected source paths: "+repr(sorted(names^SOURCES)))
        for m in members:
            path=ROOT/m.name
            path.parent.mkdir(parents=True,exist_ok=True)
            stream=tar.extractfile(m)
            if stream is None: raise RuntimeError("Missing file: "+m.name)
            path.write_bytes(stream.read())
    lakefile=ROOT/"lakefile.toml"
    txt=lakefile.read_text()
    unused='[[lean_lib]]\nname = "RequestProject"\nglobs = ["RequestProject.+"]\n'
    if unused not in txt:
        raise RuntimeError("Unexpected Lake template: manual review required")
    lakefile.write_text(txt.replace(unused,''))
    print("Imported",len(SOURCES),"Lean/toolchain/manifest files from pinned source")
    print("Excluded RequestProject template library; left existing Python sources intact")

if __name__=="__main__":
    main()
