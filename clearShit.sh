rm -rf gens_*
rm -rf gensidm_*
cd runs && \
find . -name scanner_plugins -type d -exec rm -rf {} \; && \
find . -name logs -type d -exec rm -rf {} \; && \
cd ..
