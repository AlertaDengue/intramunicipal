# Correção da geometria do Bongi

A referência original do Recife contém uma geometria inválida para o bairro Bongi (`id_bairro = 809`) em `rf_ref` e `rf_rpa_shape`. O script `fix_bongi_geometry.R` usa `sf::st_make_valid()` somente nessa geometria e preserva as demais geometrias e os atributos da referência.

É necessário ter R e o pacote `sf` instalados. Execute a partir da raiz do repositório, usando como entrada uma cópia **original e ainda não corrigida** da referência e caminhos distintos para preservar o arquivo de origem:

```sh
Rscript recife/fix_bongi_geometry.R \
  recife/ref_RF_original.RData \
  recife/ref_RF_corrigido.RData
```

O script valida o resultado antes de salvá-lo. O `recife/ref_RF.RData` corrigido já está incluído na PR #5; não é necessário executar o script durante a geração normal do boletim.
