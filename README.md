Este repositório tem como objetivo a produção dos produtos do projeto já realizados para municípios,
estados, etc., para escala espacial mais fina (intramunicipal). 

Os produtos disponíveis no boletim intramunicipal são:
1. curva de incidência com nowscast e níveis de alerta para o município;
2. curva de incidência com nowscast e níveis de alerta em escala intramunicipal (joinville -  distritos, recife - rpa);
3. comparação histórica município;
4. Rt município;
5. Rt intramunicipal;
6. Mapa de incidência intramunicipal;
7. Receptividade climática

Na raiz do diretório está o script "funcoes_submunicipal.R" que concatena toda as funções necessárias para execução
do boletim.

Até o momento, os municípios incluídos  são: Joinville, Recife e Brasília. Cada um destes
municípios possui sua própria pasta com o arquivo .qmd que gera o boletim de cada município individual,
um RData com os arquivos de referência necessários (código de identificação dos bairros,
shapefiles da escala intramunicipal, população dos bairros). Durante a execuação do .qmd,
o arquivo "cabeçalho_municipio.tex" para incluir informação da SE vigente no cabeçalho do boletim. 

Em cada pasta de município, está o .qmd para gerar o boletim intramunicipal.
Para construir o boletim são necessários os dados abaixo de cada município:
1.tabela.historico_alerta - dengue;
2.tabela.historico_alerta - chik;
3.tabela.notificação

Para realizar testes, estes dados são baixados da API e salvos no diretório
"dados_intramunicipal_teste/joinvile/sexx". Dessa forma são efetuados eventuais debugs no código.

Também está na raiz do diretório o arquivo "boletim_intramunicipal_consolidado.qmd" que é uma proposta de
.qmd generalizável para qualquer município que entre na análise intramunicipal. A partir da identificação
do código do município, os produtos específicos para o município são gerados.


