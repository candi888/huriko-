#!/bin/sh

gfortran -Wall Calc.f90 && ./a.out
echo "fortran 計算終了"
python -u "anim.py"
echo "python 描画終了"