# Pipeline de Anotación de Variantes y Conservación de Proteínas

## Descripción General

Este pipeline implementa un flujo de trabajo reproducible para analizar la variación genética dentro de secuencias codificantes (CDS) y evaluar la conservación a nivel de proteína. El pipeline:

- Extrae regiones codificantes para genes de interés
- Identifica variantes genéticas dentro de estas regiones
- Anota los efectos a nivel de codón y aminoácido
- Reconstruye la variación proteica a partir de genotipos poblacionales
- Cuantifica frecuencias de aminoácidos por posición
- Genera visualizaciones de conservación y variación proteica

El pipeline está diseñado para ser generalizable a cualquier gen y conjunto de transcritos, dados los archivos de entrada apropiados.

---

## Flujo Conceptual

```
Genoma + GFF + VCF
        ↓
Modelado de CDS por transcrito
        ↓
Extracción de secuencias (CDS + proteína)
        ↓
Extracción de variantes (región → CDS)
        ↓
Anotación de efectos (codón → aminoácido)
        ↓
Resúmenes por sitio y transcrito
        ↓
Reconstrucción basada en genotipos
        ↓
Matriz de frecuencias de aminoácidos
        ↓
Visualización (gráficas de barras apiladas)
```

---

## Estructura del Proyecto

```
spo11_variation/
├── config/
│   └── config.yaml
├── data/
│   ├── genome/
│   ├── annotation/
│   ├── hapmap/
│   └── project_inputs/
├── results/
│   └── <nombre_proyecto>/
│       ├── 01_models/
│       ├── 02_sequences/
│       ├── 03_variants/
│       ├── 04_reports/
│       └── logs/
├── scripts/
└── workflow/
    └── Snakefile
```

---

## Configuración

Todos los parámetros del pipeline se controlan mediante `config/config.yaml`. A continuación se muestra un ejemplo de configuración:

```yaml
project_name: "spo11_agpv3"

gene:
  gene_id: "GRMZM2G129913"
  gene_name: "SPO11"

transcripts:
  - "T01"
  - "T02"
  - "T03"

paths:
  genome_fa: "data/genome/B73.fa.gz"
  annotation_gff: "data/annotation/annotation.gff3.gz"
  hapmap_vcf: "data/hapmap/input.vcf.gz"
  project_inputs_dir: "data/project_inputs"

results_dir: "results/spo11_agpv3"
```

---

## Ejecución

### Prueba (dry run)

```bash
snakemake -s workflow/Snakefile -n
```

### Ejecución completa

```bash
snakemake -s workflow/Snakefile --cores 4
```

---

## Módulos del Pipeline

### 1. Validación de entradas

**Script:** `validate_inputs.py`

Valida la existencia de archivos, consistencia de nombres de cromosomas, presencia de transcritos en la anotación y disponibilidad de CDS.

### 2. Construcción de modelos CDS

**Script:** `build_cds_bed_from_gff.py`

Genera archivos BED para extracción de secuencias e intersección de variantes.

### 3. Mapeo de coordenadas CDS

**Script:** `build_blockmap.py`

Convierte coordenadas genómicas a posiciones CDS, manejando estructuras multi-exónicas y orientación de la cadena.

### 4. Extracción de secuencias CDS

**Script:** `extract_cds_fasta.py`

Produce secuencias CDS concatenadas por transcrito y bloques individuales de exones para depuración.

### 5. Traducción a proteína

**Script:** `translate_cds_to_pep.py`

Genera secuencias proteicas con validación de longitud correcta y ausencia de codones de parada internos.

### 6. Extracción de variantes regionales

**Script:** `extract_region_variants_from_vcf.py`

Filtra el VCF completo a la región genómica del gen de interés.

### 7. Subconjunto de variantes por CDS

**Script:** `subset_variants_by_bed.py`

Extrae variantes específicas para cada transcrito dentro de los límites de CDS.

### 8. Anotación de efectos de variantes

**Script:** `annotate_cds_variants.py`

Determina el codón afectado, el aminoácido resultante y el tipo de variante (sinónima, no sinónima, nonsense).

### 9. Resumen por sitio

**Script:** `summarize_cds_sites.py`

Agrupa variantes por posición CDS y codón, resumiendo los efectos observados.

### 10. Resumen por transcrito

**Script:** `summarize_transcripts.py`

Calcula el número de variantes, proporción de efectos y longitudes de CDS/proteína por transcrito.

### 11. Extracción de genotipos

**Script:** `extract_region_genotypes_from_vcf.py`

Genera un VCF que contiene genotipos para todas las muestras dentro de la región del gen.

### 12. Subconjunto de genotipos por CDS

**Script:** `subset_genotypes_by_bed.py`

Produce VCFs específicos por transcrito que contienen solo variantes en CDS.

### 13. Matriz de frecuencias de aminoácidos

**Script:** `build_aa_frequency_matrix_from_vcf.py`

Reconstruye secuencias de aminoácidos por individuo y calcula frecuencias por posición. Genera archivos en formato largo (`*.long.tsv`) y matriz (`*.matrix.tsv`).

### 14. Visualización

**Script:** `plot_aa_frequency_stacked.py`

Genera dos gráficas de barras apiladas por transcrito:

- **Composición completa:** Todos los aminoácidos (suma = 100%)
- **Variación no-referencia:** Se elimina el aminoácido dominante para resaltar la variación

---

## Archivos de Salida

### Efectos de variantes

`03_variants/effects/` — Anotaciones de efectos por transcrito

### Resúmenes

```
04_reports/
├── cds_site_summary.tsv       # Resumen de variantes por sitio
└── transcript_summary.tsv     # Estadísticas por transcrito
```

### Frecuencias de aminoácidos

`SPO11.T03.aa_freq.matrix.tsv` — Matriz de frecuencias por posición

### Figuras

```
SPO11.T03.aa_freq.stacked.nonref.final.png   # Variación no-referencia
SPO11.T03.aa_freq.stacked.full.final.png     # Composición completa de aminoácidos
```

---

## Interpretación Biológica

Este pipeline permite:

- Identificar posiciones altamente variables
- Detectar sitios con posibles efectos funcionales
- Evaluar la conservación de proteínas en poblaciones
- Comparar isoformas alternativas de un mismo gen

---

## Supuestos Importantes

- La anotación GFF es consistente con el ensamblaje del genoma
- El VCF está alineado al mismo ensamblaje de referencia
- Los transcritos de interés están definidos en `project_inputs`
- La nomenclatura sigue el formato: `{gene_id}_{transcript}` y `{gene_name}.{transcript}`

---

## Extensiones Futuras

- Anotación de dominios funcionales
- Análisis comparativos multi-gen

---

## Reproducibilidad

El pipeline es:

- **Determinista:** Produce resultados consistentes con entradas idénticas
- **Modular:** Cada paso puede ser validado independientemente
- **Configurable:** Todos los parámetros centralizados en `config.yaml`

La reproducibilidad completa se logra a partir de las siguientes entradas:

- Ensamblaje del genoma (FASTA)
- Anotación génica (GFF)
- Datos de variación poblacional (VCF)
- Definiciones de transcritos

---

## Licencia

Este proyecto está bajo la licencia MIT. Ver el archivo [LICENSE](https://www.notion.so/LICENSE) para más detalles.

---

## Contacto

Para preguntas o problemas, por favor abre un issue en el repositorio o contacta a:

Marina Mendoza Suárez

[marinams@lcg.unam.mx](mailto:marinams@lcg.unam.mx)

---