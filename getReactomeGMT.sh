wget -N "https://reactome.org/download/current/ReactomePathways.gmt.zip"
unzip -o ReactomePathways.gmt.zip
sed -i 's/ /_/g' ReactomePathways.gmt

