cut -f1 *.beagle > *.sites
cut -d '_' -f1-2 *.sites > *.sites1
cut -d '_' -f3 *.sites > *.sites2
paste *.sites1 *.sites2 | column -s $'\t' -t > *.sites
echo "$(tail -n +2 *.sites)" > *.sites
rm *.sites1
rm *.sites2