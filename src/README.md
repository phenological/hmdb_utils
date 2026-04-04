# HMDB Local API

A FastAPI wrapper around the locally compiled HMDB metabolite and 1H-NMR databases.

## Project layout

```
hmdb_api/
├── main.py               ← FastAPI application (this file)
├── hmdb_local_tools.py   ← Core search/scoring functions (copy yours here)
├── requirements.txt
└── inst/
    ├── serum_metabolites_with_spectra.csv
    └── spectral1hnmr.csv
```

## Setup

```bash
pip install -r requirements.txt
```

Copy your `hmdb_local_tools.py` into the same directory as `main.py`.

Place your compiled CSVs under `inst/` (or point to them via env vars):

```bash
export HMDB_METABOLITES_CSV=inst/serum_metabolites_with_spectra.csv
export HMDB_NMR_CSV=inst/spectral1hnmr.csv
```

## Run

```bash
uvicorn main:app --reload
```

Interactive docs at **http://localhost:8000/docs**

---

## Endpoints

### Metabolites

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/metabolites/search?q=<name>&limit=10` | Fuzzy search by name or synonym |
| `GET` | `/metabolites/{accession}` | Full record for one metabolite |
| `GET` | `/metabolites/{accession}/description` | Free-text description only |
| `GET` | `/metabolites/{accession}/synonyms` | All synonyms |
| `GET` | `/metabolites/{accession}/formula` | Chemical formula + atom counts |

### NMR

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/nmr/{accession}` | All recorded 1H-NMR signals for a metabolite |
| `POST` | `/nmr/query` | Query DB by signal pattern (see below) |

### Utilities

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/formula/parse?formula=C6H12O6` | Parse any formula string |
| `GET` | `/health` | Confirm data is loaded |

---

## NMR query body

```json
{
  "signals": [
    {
      "range": [7.07, 7.18],
      "mult": "d",
      "ppm": [7.1228, 7.136],
      "heights": [0.83, 1.0]
    },
    {
      "range": [7.2, 7.35],
      "mult": "d",
      "ppm": [7.241, 7.2542],
      "heights": [1.0, 0.83]
    }
  ],
  "limit": 6
}
```

- `range` – `[from_ppm, to_ppm]` window where the peak must exist  
- `mult`  – multiplicity code (`s`, `d`, `t`, `q`, `m`, `dd`, …) or `"*"` for wildcard  
- `ppm` / `heights` – optional multiplet shape for finer scoring (auto-normalised)  

Results are sorted by `similarity` descending.

---

## LLM tool definitions (example)

If you are wiring this API to an LLM agent, here are compact tool stubs:

```python
tools = [
    {
        "name": "search_metabolite",
        "description": "Fuzzy search for metabolites by name or synonym.",
        "parameters": {"q": "string", "limit": "int (default 10)"},
        "http": "GET /metabolites/search?q={q}&limit={limit}",
    },
    {
        "name": "get_metabolite",
        "description": "Retrieve full details for a metabolite by HMDB accession.",
        "parameters": {"accession": "string"},
        "http": "GET /metabolites/{accession}",
    },
    {
        "name": "get_nmr_signals",
        "description": "Get all 1H-NMR signals recorded for a metabolite.",
        "parameters": {"accession": "string"},
        "http": "GET /nmr/{accession}",
    },
    {
        "name": "query_nmr",
        "description": "Find metabolites matching one or more NMR signals (AND logic). Supports multiplet-shape scoring.",
        "parameters": {"signals": "list[SignalDescriptor]", "limit": "int"},
        "http": "POST /nmr/query",
    },
]
```