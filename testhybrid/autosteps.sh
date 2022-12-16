#!/bin/sh

sed "s/StepsZenith = .*/StepsZenith = 0.5/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.5	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.90/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.90	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.95/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.95	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 1.0/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "1.0	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 1.05/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "1.05	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 1.1/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "1.1	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.15/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.15	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps


sed "s/StepsZenith = .*/StepsZenith = 0.10/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.10	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.2/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.2	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.45/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.45	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.35/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.35	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.25/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.25	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps


sed "s/StepsZenith = .*/StepsZenith = 0.55/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.55	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.40/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.40	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.30/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.30	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps


sed "s/StepsZenith = .*/StepsZenith = 0.60/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.60	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps


sed "s/StepsZenith = .*/StepsZenith = 0.65/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.65	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.7/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.7	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.75/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.75	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.8/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.8	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps

sed "s/StepsZenith = .*/StepsZenith = 0.85/g" compare_both.py > .tempfile 
cat .tempfile > compare_both.py && rm .tempfile
echo -n "0.85	" >> comparesteps
python compare_both.py
cat hybrid_sigmas >> comparesteps


