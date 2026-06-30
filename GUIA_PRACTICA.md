# GUÍA PRÁCTICA: Descargando datos de SKEMPI para tu GA

## 1. BARNASE-BARSTAR (El clásico para empezar)

### Opción A: Usar el script automático
```bash
# Descargar SKEMPI, extraer interfaz y obtener PDB
python skempi_downloader.py --pdb 1B27 --output ./barnase_barstar/

# Esto genera:
# - 1B27.pdb                              (estructura)
# - 1B27_interface_positions.txt          (posiciones de interfaz)
# - 1B27_mutations.txt                    (mutaciones experimentales)
# - 1B27_detailed_mutations.csv           (datos detallados con ddG)
```

**Output esperado:**
```
Interface residues for 1B27:
Chain A interface positions:
20 27 48 56 60 62 85 102
Chain D interface positions:
28 29 30 35 39 42 43
```

### Opción B: Cálculo manual de distancia
```bash
python extract_interface_residues.py \
  -pdb 1B27.pdb \
  -c1 A \
  -c2 D \
  -method distance \
  -cutoff 5.0
```

---

## 2. ANTI-LISOZIMA D44.1 (Para antibodies)

```bash
python skempi_downloader.py --pdb 1B4F --output ./d44_lysozyme/

# Esto te dará:
# - 1B4F.pdb con >18 mutaciones experimentales
# - Posiciones de CDR donde se han hecho las mejoras de afinidad (140x)
```

---

## 3. VER TODOS LOS COMPLEJOS DISPONIBLES

```bash
python skempi_downloader.py --list-all | head -20
```

Ejemplo de salida:
```
Total de complejos en SKEMPI v2.0: 81

1A2K: 8 mutations
1ACB: 14 mutations
1AVX: 6 mutations
1B27: 65 mutations        ← BARNASE-BARSTAR (RECOMENDADO)
1B3S: 12 mutations
1B4F: 18 mutations        ← D44.1 ANTIBODY
...
```

---

## 4. PARA TU ALGORITMO GENÉTICO

Una vez que hayas descargado, puedes usar directamente:

```python
# En tu código de GA:
import pandas as pd

# Cargar las posiciones de interfaz
with open('barnase_barstar/1B27_interface_positions.txt', 'r') as f:
    lines = f.readlines()
    # Parsear Chain A positions
    chain_a_pos = list(map(int, lines[2].split()))
    # Parsear Chain D positions  
    chain_d_pos = list(map(int, lines[6].split()))

# Cargar datos experimentales de ddG
exp_data = pd.read_csv('barnase_barstar/1B27_detailed_mutations.csv')

# Tu GA puede entonces:
# 1. Mutar solo residuos en chain_a_pos
# 2. Comparar sus predicciones con exp_data['ddG']
# 3. Calcular correlación Pearson para benchmarking
```

---

## 5. ESTRUCTURA RECOMENDADA DE CARPETAS

```
tu_proyecto_GA/
├── data/
│   ├── barnase_barstar/
│   │   ├── 1B27.pdb
│   │   ├── 1B27_interface_positions.txt
│   │   ├── 1B27_mutations.txt
│   │   └── 1B27_detailed_mutations.csv
│   ├── d44_lysozyme/
│   │   ├── 1B4F.pdb
│   │   └── ...
│   └── skempi_v2.csv (caché descargado)
├── src/
│   ├── genetic_algorithm.py
│   ├── rosetta_interface.py
│   └── benchmark.py
└── results/
    ├── barnase_predictions.csv
    └── correlations.txt
```

---

## 6. TABLA RÁPIDA DE COMPLEJOS RECOMENDADOS POR CASO DE USO

### Para empezar (pequeños, rápidos)
| PDB  | Proteína        | # Mutaciones | Cadenas | Complejidad |
|------|-----------------|-------------|---------|-------------|
| 1B27 | Barnase-Barstar | 65          | A-D     | Baja        |
| 1B3S | Barnase-Barstar | 12          | A-D     | Muy baja    |
| 1G4Y | Ubiquitin-UBC   | 20          | A-B     | Media       |

### Intermedio (antibodies)
| PDB  | Proteína          | # Mutaciones | Cadenas | Complejidad |
|------|-------------------|-------------|---------|-------------|
| 1B4F | D44-Lysozyme      | 18          | H-L     | Media       |
| 1A2K | Antibody-Lysozyme | 14          | H-L     | Media       |

### Advanced (complejos enzimáticos)
| PDB  | Proteína          | # Mutaciones | Cadenas | Complejidad |
|------|-------------------|-------------|---------|-------------|
| 2HQQ | Proteasome        | 42          | A-B     | Alta        |
| 1HVH | HIV protease      | 35          | A-B     | Alta        |

---

## 7. WORKFLOW COMPLETO PARA VALIDACIÓN

```bash
#!/bin/bash

# 1. Descargar Barnase-Barstar
python skempi_downloader.py --pdb 1B27 --output ./test_case/

# 2. Tu GA predice ddG para cada mutación en 1B27_detailed_mutations.csv
python genetic_algorithm.py \
  --pdb test_case/1B27.pdb \
  --interface test_case/1B27_interface_positions.txt \
  --output predictions_1b27.csv

# 3. Comparar con experimental
python benchmark.py \
  --predictions predictions_1b27.csv \
  --experimental test_case/1B27_detailed_mutations.csv \
  --output benchmark_results.txt

# 4. Ver correlación
# Pearson > 0.6 = buena predicción ✓
# Pearson > 0.75 = excelente ✓✓
```

---

## 8. SOLUCIÓN RÁPIDA SI NO QUIERES PROGRAMAR

Si quieres usar SKEMPI directamente sin scripts:

1. Ve a https://life.bsc.es/pid/skempi2/
2. Descarga el CSV
3. Filtra por PDB que te interese (1B27 para Barnase)
4. En la columna "Location" busca residuos con "COR" o "RIM" (core/rim interface)
5. Esas son tus posiciones de interfaz

---

## 9. REFERENCIA RÁPIDA: COLUMNS EN SKEMPI CSV

```
PDB              → ID estructura (1B27)
Mutation(s)...   → Formato: A123B (cadena_posición_aa_mutante)
Affinity_Kd_wt   → Afinidad wild-type (nM)
Affinity_Kd_mut  → Afinidad mutante (nM)
DDG(kcal/mol)    → Cambio energía libre (≤0 = mejora)
Location         → COR/RIM/SUP/INT/SUR
```

---

Ahora ya tienes TODO para empezar. ¿Necesitas ayuda con algún paso específico o quieres que te prepare un template del GA para usar con estos datos?
