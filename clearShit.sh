rm -rf gens_*
cd runs
find . -name logs -type d -exec rm -rf {} \;
find . -name scanner_plugins -type d -exec rm -rf {} \;
cd ..
