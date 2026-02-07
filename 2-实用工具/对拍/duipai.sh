#!/bin/bash
# Usage: ./duipai.sh [my] [gen] [std] [spj]
# Defaults: my.cpp gen.cpp std.cpp
#
# Arguments:
#   my:  待测程序 (default: my.cpp)
#   gen: 数据生成器 (default: gen.cpp)
#   std: 标程/校验器 (default: std.cpp)
#   spj: 开启 SPJ 模式 (optional: "spj")
#
# SPJ Protocol:
#   Cmd: ./std <input> <user_out>
#   Ret: 0 = AC, !=0 = WA

M=${1:-my.cpp}
G=${2:-gen.cpp}
S=${3:-std.cpp}

# 快速编译
g++ $M -o m -O2 || exit 1
g++ $G -o g -O2 || exit 1
g++ $S -o s -O2 || exit 1

t=0
while true; do
	((t++))
	printf "\rTest: %d" $t

	./g > data.in
	./m < data.in > data.out

	if [ "$4" == "spj" ]; then
		./s data.in data.out > data.ans || break
	else
		./s < data.in > data.ans
		diff -w data.out data.ans || break
	fi
done
echo -e "\nWA!"
