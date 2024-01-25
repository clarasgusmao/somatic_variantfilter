# diretorio
cd my_samples/

## variaveis de ambiente
CHAIN="hg19ToHg38.over.chain"
HPO="$1"
SAMPLE="";

## vep

##variavel global da pasta onde se encontra o arquivo de saida
output_folder=vep_output

#rm -rf vep_output/*.filter*
mkdir -p $output_folder
chmod 777 $output_folder

#---------------------------------------------------#
# Aqui tinha a etapa do VEP                         #
# Ver: https://github.com/renatopuga/lmabrasil-hg38 #
#---------------------------------------------------#

## filter_vep
##Laço for para pegar todos os arquivos (*.vcf) que terminem em .vcf;
for filename in *.vcf; do 
  ## dentro do laço for pego o basename do arquivo.
  base_name=$(basename ${filename})
  ## expressão regular para pegar apenas valores entre underscore
  SAMPLE=$(echo "$base_name" | grep -o '_\([^_]*\)_' | sed 's/_//g')
  ## echo para ver se peguei o basename do arquivo certo e deixando explicito que iniciei a leitura do arquivo XXXX.
  echo "STARTING SAMPLE N:$SAMPLE"
  ## removido o /vep_output para segregar o arquivo de entrada do de saida.
  udocker --allow-root run --rm  -v `pwd`:`pwd` -w `pwd` ensemblorg/ensembl-vep filter_vep \
  -i liftOver_$SAMPLE\_$(basename $CHAIN .over.chain).vep.vcf \
  -filter "(MAX_AF <= 0.01 or not MAX_AF) and \
  (FILTER = PASS or not FILTER matches strand_bias,weak_evidence) and \
  (SOMATIC matches 1 or (not SOMATIC and CLIN_SIG matches pathogenic)) and (not CLIN_SIG matches benign) and \
  (not IMPACT matches LOW) and \
  (Symbol in hpo/$HPO)"  \
  --force_overwrite \
  -o $output_folder/liftOver_$SAMPLE\_$(basename $CHAIN .over.chain).vep.filter.vcf

  ## bcftools +split-vep

  # 1. criar o cabeçalho
  bcftools +split-vep -l $output_folder/liftOver_$SAMPLE\_$(basename $CHAIN .over.chain).vep.filter.vcf | \
  cut -f2  | \
  tr '\n\r' '\t' | \
  awk '{print("CHROM\tPOS\tREF\tALT\t"$0"FILTER\tTumorID\tGT\tDP\tAD\tAF\tNormalID\tNGT\tNDP\tNAD\tNAF")}' > vep_output/liftOver_$SAMPLE\_$(basename $CHAIN .over.chain).vep.filter.tsv

  # 2. adicionar as variantes
  bcftools +split-vep \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t%FILTER\t[%SAMPLE\t%GT\t%DP\t%AD\t%AF\t]\n' \
  -i 'FMT/DP>=20 && FMT/AF>=0.1' -d -A tab $output_folder/liftOver_$SAMPLE\_$(basename $CHAIN .over.chain).vep.filter.vcf \
  -p x  >> $output_folder/liftOver_$SAMPLE\_$(basename $CHAIN .over.chain).vep.filter.tsv
  echo "$SAMPLE - DONE!"
done;