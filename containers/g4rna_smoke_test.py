"""Build-time smoke test for containers/Dockerfile.g4rna.

Fails the image build if the tool's own bundled sample no longer
reproduces its own documented expectations under whichever numpy/pandas/
PyBrain versions actually resolved -- a pinned-version drift here would
otherwise surface as a silently wrong score much later, at first real
Atlas use.
"""
import subprocess
import sys

out = subprocess.check_output([
    "python", "screen.py", "sample.fas", "-a", "G4RNA_2016-11-07.pkl",
    "-c", "description", "G4NN", "-w", "60", "-s", "10", "-e",
])
rows = [line.split("\t") for line in out.splitlines()[1:]]
terra = next(r for r in rows if "TERRA" in r[1])
spinach = [r for r in rows if "false negative" in r[1]]

assert float(terra[2]) > 0.9, "TERRA (a confirmed G4) scored low: " + terra[2]
assert all(float(r[2]) < 0.3 for r in spinach), "a documented non-G4 scored high"
print("smoke test passed: %d rows, TERRA G4NN=%.3f" % (len(rows), float(terra[2])))
sys.exit(0)
