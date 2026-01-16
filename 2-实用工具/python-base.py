# 快速读入（适合大数据量）
import sys
input = lambda: sys.stdin.readline().strip()

# 高精度整数（Python内置int支持任意精度）
a = int(input())
b = int(input())
print(a + b)  # 加法
print(a * b)  # 乘法
print(a // b) # 整除
print(a % b)  # 取模

# 数组定义（列表）
arr = [int(x) for x in input().split()]
print(arr)

# 全局变量定义
global_var = 12345678901234567890

# 函数定义
def add(x, y):
	return x + y
print(add(a, b))

# 读入多行（如n行数据）
n = int(input())
for _ in range(n):
	x = int(input())
	print(x * 2)

# 常用操作补充
# 排序
arr.sort()
print(arr)

# 最大最小值
print(max(arr), min(arr))

# 字符串处理
s = input()
print(s[::-1])  # 反转字符串
print(s.upper(), s.lower())  # 大小写转换

# 字典（哈希表）
counter = {}
for x in arr:
    counter[x] = counter.get(x, 0) + 1
print(counter)

# 集合（去重）
unique = set(arr)
print(unique)

# 二维数组
matrix = [[0]*3 for _ in range(2)]
print(matrix)

# 枚举（遍历下标和元素）
for idx, val in enumerate(arr):
    print(idx, val)

# lambda表达式
f = lambda x: x**2
print(f(5))