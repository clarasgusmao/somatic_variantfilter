# Filtro e Análise de Variantes de Mielofibrose

Hi! I'm your first Markdown file in **StackEdit**. If you want to learn about StackEdit, you can read me. If you want to play with Markdown, you can edit me. Once you have finished with me, you can create new files by opening the **file explorer** on the left corner of the navigation bar.


# Geral

Este repositório contem todos os aqruivos necessários para a filtragem e análise de variantes de Mielofibrose. Nele, estão contidos os dados anotados para análise, o script de filtragem e arquivos extras.

Os códigos referenciados neste repositório foram feitos para serem utilizados em notebooks jupyter, ou google colab.

## 1. Clonando o Repositório 

Para clonar este repositório utilize o código:

```
%%bash
https://github.com/clarasgusmao/somatic_variantfilter.git
```

## 2. Instalação das Ferramentas Necessárias

### a) Instalação do BCFtools com plugin split-vep
O plugin permite extrair os campos de anotações estruturadas como INFO/CSQ criadas por bcftools/csq ou VEP (em nosso caso VEP).

Mais informações: https://samtools.github.io/bcftools/howtos/plugin.split-vep.html
```
%%bash
git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
make install
```

### b) Instalação do uDocker

Udocker é uma ferramenta básica para executar containers docker simples em sistemas sem privilégios de root. Esta abordagem de rodar como sudo será apenas utilizada no workflow pois ele foi elaborado para ser utilizado com o Google Colab e não é recomendada.

Em nosso caso, sempre que utilizarmos o comando udocker rodamos com a opção  `docker --allow-root`.

Mais informações: [https://indigo-dc.github.io/udocker/](https://indigo-dc.github.io/udocker/)

```
%%bash
pip install udocker
udocker --allow-root install
```

### c) Download da imagem do ensembl-vep
Ensembl vep é um conjunto de ferramentas para predição de impactos de variantes. Neste workflow, usaremos o comando de filtragem do vep, para filtrar as variáveis de interesse. Como o udocker foi instalado, é possível fazer o download da imagem do vep usando udocker --allow-root pull.

Mais informações: https://grch37.ensembl.org/info/docs/tools/vep/index.html

```
%%bash
udocker --allow-root pull ensemblorg/ensembl-vep
```


## 3. Filtragens das Variantes de Interesse

A filtragem das variáveis é executada no script vep-gc.sh. O arquivo Myelofibrosis.txt possui uma lista dos genes de interesse. Para mais informações sobre o script, favor consultar o script completo no repositório.

(Neste workflow, foram utilizadas as amostras do projeto LMA Brasil. Os arquivos VCF do projeto foram convertidos previamente da versão do genoma hg19 para hg38 utilizando o programa gatk LiftoverVcf com as posições hg19ToHg38.over.chain da UCSC.)

```
%%bash
sh somatico_2024/vep-gc.sh Myelofibrosis.txt
```
- Filtragem de variantes para definição de Score MSDmanual
```
%%bash
sh somatico_2024/vep-gc.sh msdgenes.txt
```

- Filtragem de variantes importantes mencionadas em artigos de revisão
```
%%bash
sh somatico_2024/vep-gc.sh revisaogenes.txt
```

## 4. Análise dos Resultados

### a) Gerar uma tabela unificada com as variáveis de interesse a partir das saídas do script vep-gc.sh. 

Mais informações:  [https://pandas.pydata.org/](https://pandas.pydata.org/)

```
import pandas as pd
import os

pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', None)

output_folder = '/content/somatico_2024/my_samples/vep_output/'
output_path = '/content/df_merged.csv'

dataframes = []

for final_outputs in os.listdir(output_folder):
    if final_outputs.endswith('.tsv'):

        output_pathway = os.path.join(output_folder, final_outputs)
        df = pd.read_csv(output_pathway, sep='\t', index_col=False)
        df['Sample'] = final_outputs.split('_')[1]
        dataframes.append(df)

df_merged = pd.concat(dataframes, ignore_index=True)
df_merged.to_csv(output_path, index=False)

```

### b) Gráfico que mostra o número de variantes que cada amostra apresenta. 
As amostras que não aparecem no gráfico não possuem variantes de interesse.

```
df_merged_sample_chart = df_merged.value_counts("Sample")
df_merged_sample_chart.plot.bar(x='Sample')
```

### c) Gerar gráfico de distribuição das variantes de interesse para cada gene. 
Alguns genes não aparecem no gráfico pois não foram encontradas variantes de interesse nos mesmos, dentre as amostras analisadas.
```
df_merged_chr_chart = df_merged.value_counts("SYMBOL")
df_merged_chr_chart.plot.pie(y='SYMBOL', figsize=(5, 5), autopct='%1.1f%%', startangle=90)
```

### d) Gráfico que relaciona a patogenicidade das variantes de acordo com os genes estudados. 
Alguns genes não aparecem no gráfico pois não foram encontradas variantes de interesse nos mesmos, dentre as amostras analisadas.

Dos genes em questão, encontramos anotações de patogenicidade para os genes KRAS, NRAS e TP53.
```
clin_sig_split = pd.DataFrame()
clin_sig_split_and = pd.DataFrame()
clin_sig_split_or = pd.DataFrame()

clin_sig_split_and = df_merged['CLIN_SIG'].str.split('&',expand=True)
for col1 in clin_sig_split_and:
  if clin_sig_split_and[col1].str.contains('/').any():
    clin_sig_split_or = clin_sig_split_and[col1].str.split('/',expand=True)
    for col2 in clin_sig_split_or:
      clin_sig_split['CLIN_SIG_'+str(len(clin_sig_split_and.columns)+col2+1)] = clin_sig_split_or[col2]
  else:
    clin_sig_split['CLIN_SIG_'+str(col1)] = clin_sig_split_and[col1]

clin_sig_split_pathogenic = clin_sig_split.where(clin_sig_split=='pathogenic')
clin_sig_split_likely_pathogenic = clin_sig_split.where(clin_sig_split=='likely_pathogenic')

df_merged_clin_chart = pd.DataFrame()
df_merged_clin_chart['SYMBOL'] = df_merged['SYMBOL']
df_merged_clin_chart['Pathogenic'] = clin_sig_split_pathogenic[clin_sig_split_pathogenic.columns[0:]].apply(lambda x: int(bool(''.join(x.dropna().astype(str)))),axis=1)
df_merged_clin_chart['Likely Pathogenic'] = clin_sig_split_likely_pathogenic[clin_sig_split_likely_pathogenic.columns[0:]].apply(lambda x: int(bool(''.join(x.dropna().astype(str)))),axis=1)
df_merged_clin_chart['Likely Pathogenic'] = df_merged_clin_chart['Likely Pathogenic'] - ( df_merged_clin_chart['Pathogenic'] * df_merged_clin_chart['Likely Pathogenic'] )

df_merged_clin_chart[( df_merged_clin_chart['Pathogenic'] != 0 ) | ( df_merged_clin_chart['Likely Pathogenic'] != 0 )].groupby('SYMBOL').sum().plot.bar(stacked=True, title='Pathogenicity of variants per gene')
```


# Conclusão

xxx
