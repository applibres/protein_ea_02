# Resumen del proyecto

## Vision general

Este repositorio implementa algoritmos evolutivos para diseno de proteinas, con foco en mutar una cadena/interfaz de una estructura PDB y evaluar variantes mediante Rosetta/PyRosetta. El proyecto combina:

- Rosetta y PyRosetta para mutacion estructural, relajacion local y calculo de metricas de interfaz.
- ESM2, via `transformers`, para proponer reemplazos de aminoacidos y evaluar log-likelihood de secuencias.
- DEAP para el modo evolutivo mono-objetivo.
- Un flujo MOEA/D propio, apoyado en `pymoo`, para optimizacion multiobjetivo.
- ProteinMPNN incluido en el repositorio como herramienta/modelo externo de diseno de secuencias.

El punto de entrada principal es `eaprot_call.py`, normalmente invocado desde `protmut_run.sh`.

## Estructura principal

- `README.md`: instrucciones basicas de instalacion, configuracion y uso.
- `Dockerfile`: construye una imagen sobre `rosettacommons/rosetta`, instala dependencias Python y deja `/app` como directorio de trabajo.
- `requirements.txt`: dependencias Python principales: `numpy`, `deap`, `pandas`, `matplotlib`, `torch`, `transformers`, `pymoo`, `Bio`.
- `eaprot_call.py`: parsea argumentos CLI, crea directorio de salida y despacha al algoritmo mono-objetivo o multiobjetivo.
- `protmut_run.sh`: wrapper shell para ejecutar `eaprot_call.py`.
- `replicates.sh`: ejecuta multiples replicas variando la semilla.
- `prot_GA.py`: algoritmo genetico mono-objetivo basado en DEAP.
- `prot_mob_pymoo.py`: algoritmo multiobjetivo tipo MOEA/D con mutacion, crossover, seleccion por vecindario y frente de Pareto.
- `utils.py`: funciones auxiliares para guardar poblaciones, estadisticas, archivos `.sc`, arbol genealogico JSON y distancias de Hamming.
- `prot_interface/`: capa de integracion con Rosetta, PyRosetta, ESM2 y configuracion del problema.
- `genetic_operators/`: representacion de individuos y operadores MOEA/D.
- `scenarios/`: casos de prueba con PDBs, archivos de caras/interfaz y matrices MSA.
- `ProteinMPNN/`: copia local de ProteinMPNN, pesos, ejemplos, notebooks y scripts auxiliares.

## Flujo de ejecucion

1. El usuario llama `protmut_run.sh` con escenario, parametros de simulacion, parametros evolutivos, indices de fitness, checkpoint, modo multiobjetivo y semilla.
2. `eaprot_call.py` valida los 9 argumentos esperados, parsea strings tipo `key=value`, prepara `/output/run<seed>` y decide el algoritmo:
   - `mobj=False`: usa `prot_GA.deap_sga_protein`.
   - `mobj=True`: usa `prot_mob_pymoo.pymoo_sga_protein`.
3. El algoritmo instancia `prot_interface.prot_problemI.prot_problem`, que funciona como fachada del problema biologico.
4. `prot_problem` extrae posiciones relevantes de la cadena ligando/interfaz, muta estructuras, evalua fitness y puede aplicar crossover estructural.
5. Las variantes se guardan como PDBs por generacion (`g0`, `g1`, etc.), junto con scorefiles `.sc`, representaciones de secuencia, CSVs y checkpoints.

Ejemplo documentado en el README:

```bash
./protmut_run.sh test04 sea \
  pdbfile=protein01.pdb,partners=A_C,ligand_chain=C \
  gen=5,popsize=5,mutp=1 \
  fitness_idxs=2,3 \
  checkpoint=False checks=2 mobj=True randomseed=15
```

## Componentes de `prot_interface`

### `prot_settingsI.py`

Define rutas y parametros globales:

- `ROSETTA_BIN`: directorio de binarios Rosetta.
- `CONFIG_PATH`: ruta base de escenarios, actualmente `/app/scenarios/`.
- `INTERF_EN`: binario Rosetta para energia de interfaz.
- `INTERF_AN`: binario Rosetta `InterfaceAnalyzer`.
- `FACE1_FILE_NAME`, `FACE2_FILE_NAME`: archivos que describen las caras de interfaz.
- `MSA_MATRIX`: nombre esperado de matriz MSA.
- `SCORE_INDEXES`: indices extraidos del scorefile de Rosetta.
- `SCORE_OBJECTIVE`: direccion de optimizacion de cada score (`1` maximizar, `-1` minimizar).
- `FLAGS`: flags usados con `InterfaceAnalyzer`.

### `prot_problemI.py`

Es la clase central del problema:

- Construye objetos de energia (`prot_energyInterf`), mutacion (`prot_mut`) y extraccion de aminoacidos (`prot_aa_extract`).
- Genera archivos de energia de interaccion y mapea posiciones relativas/absolutas de residuos.
- Crea el individuo inicial desde el PDB del escenario.
- Extrae secuencias completas o solo posiciones de interes.
- Evalua fitness individual o poblacional en paralelo con `multiprocessing`.
- Muta poblaciones en paralelo.
- Aplica crossover estructural en paralelo.

Las metricas comentadas para fitness incluyen `packstat`, `sc_value`, `total_score`, `delta_unsatHbonds`, `fa_rep`, `per_residue_energy_int`, `dSASA_int`, `dG_separated/dSASAx100` y `hbonds_int`.

### `prot_mutI.py`

Implementa mutacion estructural:

- Inicializa PyRosetta.
- Lee posiciones mutables desde el archivo `FACE2_FILE_NAME`.
- Ejecuta `rosetta_scripts.default.linuxgccrelease` con `PM_Mutation_Relax_Local.xml`.
- Selecciona posiciones a mutar segun energias de interaccion.
- Usa ESM2 para elegir reemplazos plausibles, con temperatura decreciente por generacion.
- Encadena mutaciones temporales y mueve el PDB final al destino.
- Implementa `crossover`, aplicando una lista de cambios posicion-aminoacido sobre una estructura base.

### `prot_aa_stI.py`

Extrae aminoacidos estabilizados/no estabilizados:

- Ejecuta `interface_energy`.
- Extrae el bloque `##### PAIRWISE SHORT-RANGE ENERGIES #####`.
- Filtra interacciones donde participa la cadena ligando.
- Agrupa energias por residuo.
- Clasifica residuos con energia total positiva como no estabilizados y negativa como estabilizados.

### `prot_energyInterfI.py`

Evalua una estructura con Rosetta `InterfaceAnalyzer`, lee el scorefile `.sc` y extrae las columnas definidas en `SCORE_INDEXES`.

### `prot_esm2.py`

Wrapper de ESM2 (`facebook/esm2_t6_8M_UR50D`):

- Selecciona dispositivo `cuda`, `mps` o CPU.
- Calcula log-likelihood total de una secuencia.
- Genera matriz de probabilidades por posicion.
- Propone reemplazos excluyendo el aminoacido wild-type y usando sampling con temperatura.
- Cachea forward passes por secuencia.

### `prot_pyrosettaI.py`

Contiene utilidades PyRosetta antiguas o complementarias para secuencia, mutacion puntual, packing, unbinding y calculo de energia bound/unbound.

## Algoritmo mono-objetivo (`prot_GA.py`)

Usa DEAP:

- Crea individuos como listas de aminoacidos con metadatos (`id`, `father`, `pdb`, `nmut`).
- Evalua un unico indice de fitness.
- Genera poblacion inicial mutando el PDB original.
- Usa elitismo del 10% y seleccion por torneo.
- Muta toda la poblacion en cada generacion.
- Guarda PDBs, scorefiles, poblacion, estadisticas y checkpoint.

Nota: en el archivo hay referencias inconsistentes a `self.fitness_idxs` aunque el constructor guarda `self.fitness_idx`. Conviene revisar este punto antes de usar intensivamente el modo mono-objetivo.

## Algoritmo multiobjetivo (`prot_mob_pymoo.py`)

Implementa un flujo MOEA/D:

- El numero de objetivos es `len(fitness_idxs) + 1`; el objetivo adicional es el delta de log-likelihood de ESM2 respecto al padre/original.
- Inicializa ESM2 al construir el algoritmo.
- Crea poblacion inicial desde el PDB original y mutaciones del individuo inicial.
- Evalua Rosetta + ESM2.
- Usa direcciones de referencia Das-Dennis, vecindarios y descomposicion Tchebycheff.
- Aplica crossover uniforme entre secuencias de vecinos y luego mutacion.
- Actualiza supervivencia por vecindario MOEA/D.
- Mantiene un frente de Pareto y lo guarda como Hall of Fame.
- Guarda checkpoint, poblacion por generacion, secuencias, scorefiles agregados y estadisticas multiobjetivo.

## Operadores geneticos

- `genetic_operators/individual.py`: representa secuencias como arrays numericos de aminoacidos (`X`), junto con fitness, PDB, padre y numero de mutaciones.
- `genetic_operators/moead.py`: contiene:
  - `CrossoverOperatorMOEAD`: crossover uniforme discreto entre dos secuencias.
  - `SelectionOperatorMOEAD`: direcciones de referencia, vecindarios y reemplazo por Tchebycheff.
- `genetic_operators/crossoever.py`: funcion experimental que usa Ollama para proponer una secuencia hija a partir de dos padres y sus scores. El nombre del archivo parece tener un typo (`crossoever`).

## Entradas esperadas por escenario

Cada escenario suele contener:

- Un PDB base.
- Archivos de caras/interfaz (`faceA.txt`, `faceB.txt`, `faceC.txt` o similares).
- Archivos `.txt` y `.sc` derivados o preexistentes.
- Opcionalmente matrices MSA.

Escenarios presentes:

- `test04`: `protein01.pdb`, `faceA.txt`, `faceC.txt`, `MSA_matrix.tsv`.
- `test05`: `6M0J_nowaters.pdb`, `faceA.txt`, `faceE.txt`, matrices MSA y mutantes de ejemplo.
- `test06`: `9Q1V_prepared_clean_relaxed.pdb`, `faceA.txt`, `faceB.txt`.

El valor por defecto de `FACE2_FILE_NAME` en configuracion es `faceB.txt`, por lo que escenarios como `test04` o `test05` pueden requerir ajustar `prot_settingsI.py` antes de correr.

## Salidas generadas

El flujo escribe principalmente en `/output/run<randomseed>`:

- `g0/`, `g1/`, ...: PDBs por generacion.
- `tmp/`: archivos temporales de mutacion/crossover.
- `checkpoint.pkl`: estado para reanudar.
- `pop_g<N>_AA.txt`: secuencias y fitness por generacion.

Tambien usa `/output/statistics` para:

- `run<seed>_individuals.csv`: poblacion historica.
- `run<seed>_individuals.json`: arbol genealogico.
- `run<seed>_scfile.csv`: agregacion de scorefiles Rosetta.
- `run<seed>_statistics.csv`: estadisticas por generacion.
- `run<seed>_hallofame.csv`: frente de Pareto en modo multiobjetivo.

El logging se escribe en `loggest.txt` y tambien en consola.

## Dependencias externas importantes

- Imagen/base Rosetta: `rosettacommons/rosetta`.
- Binarios Rosetta esperados: `InterfaceAnalyzer.default.linuxgccrelease`, `interface_energy.default.linuxgccrelease`, `rosetta_scripts.default.linuxgccrelease`.
- Base de datos Rosetta en `/usr/local/database`.
- PyRosetta importable desde Python.
- Modelo Hugging Face `facebook/esm2_t6_8M_UR50D`; si no esta cacheado, requiere acceso a red para descargar tokenizer/pesos.
- GPU opcional para acelerar ESM2.

## Observaciones y posibles riesgos

- `git` no esta instalado en el entorno actual, por lo que no se pudo inspeccionar estado de cambios con `git status`.
- El sandbox del entorno no pudo ejecutar comandos sin permisos elevados por una limitacion de namespaces; la inspeccion se hizo con comandos de lectura aprobados.
- `prot_GA.py` parece tener inconsistencias: usa `self.fitness_idxs` aunque el constructor define `self.fitness_idx`; tambien llama `mutate_population` con tuplas de 3 argumentos, mientras `prot_problem.mutate` espera 5 argumentos (`pdb_file`, `output_file`, `mut_rate`, `generation`, `ngen`).
- `replicates.sh` usa `fitness_idsx` en vez de `fitness_idxs`, lo que probablemente rompe el parsing esperado.
- `protmut_run.sh` valida ` $# -lt 8`, pero realmente usa 9 argumentos.
- Hay archivos y comentarios con mezcla de ingles/espanol y algunos typos (`algoritm`, `FITNSESSS_IDXS`, `crossoever`).
- Algunas llamadas a Rosetta usan `shell=True` y concatenacion de strings; si se aceptan rutas/parametros externos, conviene endurecer esa capa.
- El modo multiobjetivo parece ser el camino mas actualizado y consistente del proyecto.

## Resumen corto para nuevos contribuidores

El proyecto toma una proteina PDB de `scenarios/`, identifica residuos relevantes de una interfaz, genera variantes por mutacion/crossover, relaja localmente cada variante con Rosetta, evalua metricas de interfaz con Rosetta y usa ESM2 para guiar mutaciones y aportar un objetivo de plausibilidad de secuencia. El modo recomendado parece ser `mobj=True`, que ejecuta `prot_mob_pymoo.py` con MOEA/D y produce poblaciones, estadisticas, scorefiles y frente de Pareto en `/output`.
