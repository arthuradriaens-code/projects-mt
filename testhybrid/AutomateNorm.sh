#!/bin/sh

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 0.001/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "0.001	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 0.01/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "0.01	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 0.1/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "0.1	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 1/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 2/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "2	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 2.5/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "2.5	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 1.5/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1.5	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 0.5/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "0.5	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

