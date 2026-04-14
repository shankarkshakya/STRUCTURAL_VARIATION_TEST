git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug
make

ls -1 *.vcf > vcf_list.txt
./SURVIVOR/Debug/SURVIVOR merge vcf_list.txt 1000 1 1 0 0 50 merged.vcf

meaning : SURVIVOR merge <file_list> <max_distance> <type> <strands> <estimate> <no_ambiguous> <min_support> <output.vcf>

1000 → max breakpoint distance (1 kb tolerance)
1 → require same SV type (DEL with DEL, etc.)
1 → require same strand orientation
0 → do not estimate distance (use coordinates)
0 → allow ambiguous matches
1 → minimum supporting files
