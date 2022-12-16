#!/bin/sh

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 0.0001/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "0.0001	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison


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
sed "s/NormalScale = .*/NormalScale = 5/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "5	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 10/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "10	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 100/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "100	0.001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison


# Back to 1
cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/NormalScale = .*/NormalScale = 1/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile

# Let's now vary ztol
sed "s/ScaleZtol = .*/ScaleZtol = 10/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.01	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison


cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 20/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.02	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison


cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 5/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.005	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison


cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 2/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.002	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 0.1/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.0001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 0.01/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.00001	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 50/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.05	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 100/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	0.1	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

cd /home/arthur/Documents/thesis/programs/NuRadioMC/NuRadioMC/SignalProp
sed "s/ScaleZtol = .*/ScaleZtol = 1000/g" radioproparaytracing.py > .tempfile 
cat .tempfile > radioproparaytracing.py && rm .tempfile
cd /home/arthur/Documents/thesis/projects-mt/testhybrid
echo -n "1	1	" >> comparison
python compare_both.py
cat hybrid_sigmas >> comparison

