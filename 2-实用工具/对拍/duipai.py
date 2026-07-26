#!/usr/bin/env python3
"""Cross-platform duipai (对拍) script — replaces duipai.sh + duipai.bat.

Usage: duipai.py [my] [gen] [std] [spj]

Arguments:
  my:   待测程序 (default: my.cpp)
  gen:  数据生成器 (default: gen.cpp)
  std:  标程/校验器  (default: std.cpp)
  spj:  SPJ 模式，s 接收 input 和 user_out 两个参数，返回 0=AC (optional: "spj")

SPJ Protocol:
  Cmd: ./s <input> <user_out>
  Ret: 0 = AC, !=0 = WA
"""

import sys
import subprocess
from pathlib import Path

EXE_SUFFIX = ".exe" if sys.platform == "win32" else ""

M = sys.argv[1] if len(sys.argv) > 1 else "my.cpp"
G = sys.argv[2] if len(sys.argv) > 2 else "gen.cpp"
S = sys.argv[3] if len(sys.argv) > 3 else "std.cpp"
SPJ = len(sys.argv) > 4 and sys.argv[4] == "spj"

# Compile
for src, out in [(M, "m"), (G, "g"), (S, "s")]:
    if subprocess.run(["g++", src, "-o", out + EXE_SUFFIX, "-O2"]).returncode:
        sys.exit(1)

t = 0
while True:
    t += 1
    print(f"\rTest: {t}", end="", flush=True)

    # gen > data.in
    subprocess.run([f"./g{EXE_SUFFIX}"], stdout=open("data.in", "w"))

    # my < data.in > data.out
    subprocess.run([f"./m{EXE_SUFFIX}"],
                   stdin=open("data.in"), stdout=open("data.out", "w"))

    if SPJ:
        # SPJ: s data.in data.out → exit code 0 = AC
        if subprocess.run([f"./s{EXE_SUFFIX}", "data.in", "data.out"]).returncode:
            break
    else:
        # std < data.in > data.ans
        subprocess.run([f"./s{EXE_SUFFIX}"],
                       stdin=open("data.in"), stdout=open("data.ans", "w"))
        # Compare (pure Python, no diff/fc dependency)
        if Path("data.out").read_text() != Path("data.ans").read_text():
            break

print("\nWA!")
