# Database schema notes

## Goal
The CSV contains longitudinal immune cell count data from a clinical-trial dataset.

Each row in `cell-count.csv` represents **one biological sample** from one subject at one timepoint.

The real structure of the data is:

- one **project** has many **subjects**
- one **subject** has many **samples**
- one **sample** has counts for many immune cell **populations**

To model that cleanly in SQLite, the database can be normalized into five tables.

---

## Tables

### 1. `projects`
Stores one row per project.

Columns:
- `project_id` - integer primary key
- `project_code` - text unique, such as `prj1`

Why:
- avoids repeating project names in many rows

---

### 2. `subjects`
Stores one row per subject within a project.

Columns:
- `subject_id` - integer primary key
- `project_id` - foreign key to `projects`
- `subject_code` - subject identifier from the CSV, such as `sbj000`
- `condition` - e.g. melanoma, carcinoma, healthy
- `age`
- `sex`
- `treatment`
- `response` - nullable because some rows have missing response values

Constraint:
- `UNIQUE(project_id, subject_code)`

Why:
- subject metadata is constant across time for each `(project, subject)` pair in this dataset
- storing it once is cleaner than repeating it for every sample

---

### 3. `samples`
Stores one row per biological sample.

Columns:
- `sample_id` - integer primary key
- `sample_code` - unique sample identifier from the CSV
- `subject_id` - foreign key to `subjects`
- `sample_type` - e.g. PBMC or WB
- `time_from_treatment_start` - e.g. 0, 7, 14

Why:
- one subject can have multiple samples across time

---

### 4. `populations`
Stores the immune cell population names.

Columns:
- `population_id` - integer primary key
- `population_name` - unique text value

Rows inserted:
- `b_cell`
- `cd8_t_cell`
- `cd4_t_cell`
- `nk_cell`
- `monocyte`

Why:
- makes the schema scalable if more populations are added later

---

### 5. `cell_counts`
Stores one row per sample-population measurement.

Columns:
- `sample_id` - foreign key to `samples`
- `population_id` - foreign key to `populations`
- `count` - integer count

Primary key:
- `(sample_id, population_id)`

Why:
- this is the core measurement table
- each sample has one count for each immune population

---

## Why this schema is a good fit

This design separates:
- subject-level metadata
- sample-level metadata
- measurement-level data

That makes the database:
- cleaner
- easier to query
- easier to scale to many projects, subjects, samples, and populations
- better for downstream analytics than storing the raw CSV as one flat table

---

## Validations performed before loading

Before inserting data, `load_data.py` checks:

1. all required columns are present
2. `sample` is unique
3. subject metadata is consistent within each `(project, subject)` pair for:
   - `condition`
   - `age`
   - `sex`
   - `treatment`
   - `response`
4. immune cell counts are numeric and non-negative
5. `time_from_treatment_start` is numeric

If any of these checks fail, the script stops with a clear error message.

---

## Output
Running:

```bash
python load_data.py
