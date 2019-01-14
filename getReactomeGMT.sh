wget "https://reactome.org/download/current/ReactomePathways.gmt.zip"
unzip ReactomePathways.gmt.zip
sed -i 's/ /_/g' ReactomePathways.gmt
sed -i 's/\t/_/2' ReactomePathways.gmt

