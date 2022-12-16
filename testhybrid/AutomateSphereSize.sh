#!/bin/sh

notify-send "Lesgoo 🎓"
sed -i "s/SphereSize = .*/SphereSize = 50./g" compare_both.py
echo "---\n" >> comparesteps
echo "50.\n" >> comparesteps
echo "---\n" >> comparesteps
./autosteps.sh
notify-send "SphereSize 50 Done 🦆"

sed -i "s/SphereSize = .*/SphereSize = 60./g" compare_both.py
echo "---\n" >> comparesteps
echo "60.\n" >> comparesteps
echo "---\n" >> comparesteps
./autosteps.sh
notify-send "SphereSize 60 Done 🐇"

sed -i "s/SphereSize = .*/SphereSize = 70./g" compare_both.py
echo "---\n" >> comparesteps
echo "70.\n" >> comparesteps
echo "---\n" >> comparesteps
./autosteps.sh
notify-send "SphereSize 70 Done 🦌"

sed -i "s/SphereSize = .*/SphereSize = 80./g" compare_both.py
echo "---\n" >> comparesteps
echo "80.\n" >> comparesteps
echo "---\n" >> comparesteps
./autosteps.sh

notify-send "SphereSize 80 Done 📡"
