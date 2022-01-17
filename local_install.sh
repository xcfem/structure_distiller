#!/bin/sh
touch installed_files_history.txt
sudo python setup.py install --prefix=/usr/local --record installed_files.txt
echo "Updating installed files history."
cat installed_files_history.txt installed_files.txt | sort | uniq > tmp.txt
mv tmp.txt installed_files_history.txt
sudo rm installed_files.txt
