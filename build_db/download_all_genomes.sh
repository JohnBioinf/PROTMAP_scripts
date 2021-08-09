function download_genome {
	name=$1
	link=$2
	genome_dir=$3
	download_file="$(awk -F '/' '{print $NF}' <<< "$link")"
	wget "$link" -P "$genome_dir"
	gunzip "$genome_dir/$download_file"
	mv "$genome_dir/${download_file%.gz}" "$genome_dir/$name.gbk"
}

genome_dir="$data_dir/genome"

if ! [ -d "$genome_dir" ]; then
       mkdir "$genome_dir"
fi

name="anaero"
link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/131/675/GCA_014131675.1_ASM1413167v1/GCA_014131675.1_ASM1413167v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="bact"
link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/131/755/GCA_014131755.1_ASM1413175v1/GCA_014131755.1_ASM1413175v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="bifi"
link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/525/GCA_000007525.1_ASM752v1/GCA_000007525.1_ASM752v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="blautia"
link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/131/715/GCA_014131715.1_ASM1413171v1/GCA_014131715.1_ASM1413171v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="clostri"
link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/131/795/GCA_014131795.1_ASM1413179v1/GCA_014131795.1_ASM1413179v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="ecoli"
link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="ery"
link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/131/695/GCF_014131695.1_ASM1413169v1/GCF_014131695.1_ASM1413169v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"

name="lacto"
link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/014/131/735/GCA_014131735.1_ASM1413173v1/GCA_014131735.1_ASM1413173v1_genomic.gbff.gz"

download_genome "$name" "$link" "$genome_dir"
