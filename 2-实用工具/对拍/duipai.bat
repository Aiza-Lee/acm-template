@echo off
:: windows 比赛环境中通常会有 git bash，如果有的话，还是sh格式的对拍方案更加简洁
:: Usage: duipai.bat [my] [gen] [std] [spj]
:: Defaults: my.cpp gen.cpp std.cpp

set M=%1
if "%M%"=="" set M=my.cpp
set G=%2
if "%G%"=="" set G=gen.cpp
set S=%3
if "%S%"=="" set S=std.cpp
set MODE=%4

echo Compiling...
g++ %M% -o m.exe -O2
if errorlevel 1 pause & exit
g++ %G% -o g.exe -O2
if errorlevel 1 pause & exit
g++ %S% -o s.exe -O2
if errorlevel 1 pause & exit

set t=0
:loop
set /a t+=1
echo Test: %t%

g.exe > data.in
m.exe < data.in > data.out

if "%MODE%"=="spj" (
	:: SPJ Mode: s.exe input user_out
	s.exe data.in data.out > data.ans
	if errorlevel 1 goto error
) else (
	:: Normal Mode: diff with std output
	s.exe < data.in > data.ans
	fc /W data.out data.ans > nul
	if errorlevel 1 goto error
)

goto loop

:error
echo WA!
pause
