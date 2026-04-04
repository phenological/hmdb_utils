"""
HMDB Local API
==============
FastAPI application exposing metabolite and 1H-NMR lookup capabilities
built on top of the HMDB compiled CSVs and hmdb_local_tools.py.

Start with:
    uvicorn main:app --reload

Data files expected (configurable via env vars):
    HMDB_METABOLITES_CSV  – path to <subset>_metabolites_with_spectra.csv
    HMDB_NMR_CSV          – path to spectral1hnmr.csv
"""

from __future__ import annotations

import ast
import os
from contextlib import asynccontextmanager
from typing import Annotated, Optional

import numpy as np
import pandas as pd
import uvicorn
from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, Field

from hmdb_local_tools import approximate_lookup, count_atoms, multiple_query

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

METABOLITES_CSV = os.getenv("HMDB_METABOLITES_CSV", "inst/serum_metabolites_with_spectra.csv")
NMR_CSV = os.getenv("HMDB_NMR_CSV", "inst/spectral1hnmr.csv")


# ---------------------------------------------------------------------------
# Application state
# ---------------------------------------------------------------------------

class AppState:
    df: pd.DataFrame
    nmrdb: pd.DataFrame


state = AppState()


def _load_data() -> None:
    """Load and pre-process both CSVs into module-level state."""

    # --- Metabolite table ---------------------------------------------------
    df = pd.read_csv(METABOLITES_CSV)
    df["synonyms"] = df["synonyms"].apply(ast.literal_eval)
    df["synonyms_cat"] = (df["name"] + " " + df["synonyms"].str.join("")).str.lower()
    state.df = df

    # --- NMR peak table -----------------------------------------------------
    nmrdb = pd.read_csv(NMR_CSV)

    def _parse_array(x: str) -> np.ndarray:
        return np.array(ast.literal_eval(x.replace("'", "")))

    def _parse_heights(x: str) -> np.ndarray:
        arr = _parse_array(x)
        total = arr.sum()
        return arr / total if total != 0 else arr

    nmrdb["ppm"] = nmrdb["ppm"].apply(_parse_array)
    nmrdb["heights"] = nmrdb["heights"].apply(_parse_heights)

    # Sort by shift so bisect-based range queries work correctly
    nmrdb = nmrdb.sort_values("shift1(ppm)").reset_index(drop=True)
    state.nmrdb = nmrdb


# ---------------------------------------------------------------------------
# Lifespan (replaces deprecated @app.on_event)
# ---------------------------------------------------------------------------

@asynccontextmanager
async def lifespan(app: FastAPI):
    _load_data()
    yield


# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = FastAPI(
    title="HMDB Local API",
    description=(
        "Search metabolites by name/synonym and query 1H-NMR peaks "
        "against the locally compiled HMDB database."
    ),
    version="1.0.0",
    lifespan=lifespan,
)


# ---------------------------------------------------------------------------
# Pydantic schemas
# ---------------------------------------------------------------------------

class MetaboliteBasic(BaseModel):
    accession: str
    name: str

    model_config = {"from_attributes": True}


class MetaboliteDetail(MetaboliteBasic):
    secondary_accessions: Optional[list[str]] = None
    synonyms: Optional[list[str]] = None
    chemical_formula: Optional[str] = None
    iupac_name: Optional[str] = None
    cas_registry_number: Optional[str] = None
    smiles: Optional[str] = None
    inchi: Optional[str] = None
    inchikey: Optional[str] = None
    monisotopic_molecular_weight: Optional[float] = None
    average_molecular_weight: Optional[float] = None
    description: Optional[str] = None


class NameMatch(BaseModel):
    accession: str
    name: str
    matched_synonym: str
    score: float


class NmrSignal(BaseModel):
    multiplet_id: str = Field(..., alias="multiplet")
    shift_ppm: float = Field(..., alias="shift1(ppm)")
    multiplicity: str = Field(..., alias="type")
    range_from: float = Field(..., alias="from")
    range_to: float = Field(..., alias="to")
    num_hydrogens: Optional[int] = Field(None, alias="hs")
    coupling_hz: Optional[list[str]] = None
    peak_ppm: Optional[list[float]] = None
    peak_heights: Optional[list[float]] = None

    model_config = {"populate_by_name": True}


class NmrQuerySignal(BaseModel):
    """A single signal constraint for the NMR query endpoint."""
    range: tuple[float, float] = Field(
        ...,
        description="Chemical shift range [from_ppm, to_ppm]",
        examples=[(0.94, 0.99)],
    )
    mult: str = Field(
        ...,
        description="Multiplicity code (s, d, t, q, m, dd, …) or '*' for wildcard",
        examples=["t"],
    )
    ppm: Optional[list[float]] = Field(
        None,
        description="Peak positions within the multiplet (for shape matching)",
    )
    heights: Optional[list[float]] = Field(
        None,
        description="Relative intensities for each peak in `ppm` (auto-normalised)",
    )


class NmrQueryRequest(BaseModel):
    signals: list[NmrQuerySignal] = Field(
        ...,
        description="One or more signal descriptors (AND logic across all signals)",
        min_length=1,
    )
    limit: int = Field(10, ge=1, le=100, description="Maximum hits to return")

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "signals": [
                        {
                            "range": [7.07, 7.18],
                            "mult": "d",
                            "ppm": [7.1228, 7.136],
                            "heights": [0.83, 1],
                        },
                        {
                            "range": [7.2, 7.35],
                            "mult": "d",
                            "ppm": [7.241, 7.2542],
                            "heights": [1, 0.83],
                        },
                    ],
                    "limit": 6,
                }
            ]
        }
    }


class NmrQueryResult(BaseModel):
    accession: str
    name: str
    similarity: float


class AtomCounts(BaseModel):
    formula: str
    atoms: dict[str, int]


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _row_to_metabolite_detail(row: pd.Series) -> MetaboliteDetail:
    def _safe_list(val):
        if isinstance(val, list):
            return val
        try:
            return ast.literal_eval(str(val))
        except Exception:
            return None

    return MetaboliteDetail(
        accession=row.get("accession"),
        name=row.get("name"),
        secondary_accessions=_safe_list(row.get("secondary_accessions")),
        synonyms=row.get("synonyms") if isinstance(row.get("synonyms"), list) else _safe_list(row.get("synonyms")),
        chemical_formula=row.get("chemical_formula"),
        iupac_name=row.get("iupac_name"),
        cas_registry_number=row.get("cas_registry_number"),
        smiles=row.get("smiles"),
        inchi=row.get("inchi"),
        inchikey=row.get("inchikey"),
        monisotopic_molecular_weight=row.get("monisotopic_molecular_weight"),
        average_molecular_weight=row.get("average_molecular_weight"),
        description=row.get("description"),
    )


# ---------------------------------------------------------------------------
# Routes – Metabolites
# ---------------------------------------------------------------------------

@app.get(
    "/metabolites/search",
    response_model=list[NameMatch],
    summary="Fuzzy metabolite search by name or synonym",
    tags=["Metabolites"],
)
def search_metabolites(
    q: Annotated[str, Query(description="Search term (name or synonym)", min_length=1)],
    limit: Annotated[int, Query(ge=1, le=50)] = 10,
):
    """
    Return up to `limit` metabolites whose name or any synonym best matches
    the query string, ranked by fuzzy similarity score.
    """
    matches = approximate_lookup(state.df, "synonyms_cat", q, limit=limit)
    if matches is None:
        return []
    results = []
    for row in matches:
        main_name, best_synonym, row_idx, score = row
        accession = state.df.loc[int(row_idx), "accession"]
        results.append(
            NameMatch(
                accession=accession,
                name=str(main_name),
                matched_synonym=str(best_synonym),
                score=float(score),
            )
        )
    return results


@app.get(
    "/metabolites/{accession}",
    response_model=MetaboliteDetail,
    summary="Get full metabolite record by HMDB accession",
    tags=["Metabolites"],
)
def get_metabolite(accession: str):
    """
    Retrieve all stored fields for a single metabolite identified by its
    primary HMDB accession (e.g. `HMDB0001925`).
    """
    mask = state.df["accession"] == accession
    if not mask.any():
        raise HTTPException(status_code=404, detail=f"Accession '{accession}' not found.")
    row = state.df[mask].iloc[0]
    return _row_to_metabolite_detail(row)


@app.get(
    "/metabolites/{accession}/description",
    response_model=dict,
    summary="Get metabolite description text",
    tags=["Metabolites"],
)
def get_metabolite_description(accession: str):
    """Return only the free-text description for a metabolite."""
    mask = state.df["accession"] == accession
    if not mask.any():
        raise HTTPException(status_code=404, detail=f"Accession '{accession}' not found.")
    description = state.df[mask].iloc[0].get("description", "")
    return {"accession": accession, "description": description}


@app.get(
    "/metabolites/{accession}/synonyms",
    response_model=dict,
    summary="Get all synonyms for a metabolite",
    tags=["Metabolites"],
)
def get_metabolite_synonyms(accession: str):
    """Return the name and full synonym list for a metabolite."""
    mask = state.df["accession"] == accession
    if not mask.any():
        raise HTTPException(status_code=404, detail=f"Accession '{accession}' not found.")
    row = state.df[mask].iloc[0]
    synonyms = row.get("synonyms")
    if not isinstance(synonyms, list):
        try:
            synonyms = ast.literal_eval(str(synonyms))
        except Exception:
            synonyms = []
    return {"accession": accession, "name": row["name"], "synonyms": synonyms}


@app.get(
    "/metabolites/{accession}/formula",
    response_model=AtomCounts,
    summary="Parse atom counts from a metabolite's chemical formula",
    tags=["Metabolites"],
)
def get_metabolite_formula(accession: str):
    """
    Return the chemical formula and a per-element atom count dictionary
    for the requested metabolite.
    """
    mask = state.df["accession"] == accession
    if not mask.any():
        raise HTTPException(status_code=404, detail=f"Accession '{accession}' not found.")
    formula = state.df[mask].iloc[0].get("chemical_formula", "")
    if not formula:
        raise HTTPException(status_code=422, detail="No chemical formula stored for this metabolite.")
    return AtomCounts(formula=formula, atoms=count_atoms(formula))


# ---------------------------------------------------------------------------
# Routes – NMR
# ---------------------------------------------------------------------------

@app.get(
    "/nmr/{accession}",
    response_model=list[NmrSignal],
    summary="Get all 1H-NMR signals for a metabolite",
    tags=["NMR"],
)
def get_nmr_signals(accession: str):
    """
    Return every multiplet recorded in the NMR database for the given
    metabolite accession, including peak positions and relative heights.
    """
    signals = state.nmrdb[state.nmrdb["accession"] == accession]
    if signals.empty:
        raise HTTPException(
            status_code=404,
            detail=f"No NMR data found for accession '{accession}'.",
        )

    result = []
    for _, sig in signals.iterrows():
        ppm_list = sig["ppm"].tolist() if isinstance(sig["ppm"], np.ndarray) else None
        heights_list = sig["heights"].tolist() if isinstance(sig["heights"], np.ndarray) else None
        j_raw = sig.get("j(hz)")
        coupling = None
        if j_raw is not None:
            try:
                coupling = ast.literal_eval(str(j_raw))
            except Exception:
                coupling = [str(j_raw)]

        result.append(
            NmrSignal(
                **{
                    "multiplet": sig["multiplet"],
                    "shift1(ppm)": float(sig["shift1(ppm)"]),
                    "type": sig["type"],
                    "from": float(sig["from"]),
                    "to": float(sig["to"]),
                    "hs": int(sig["hs"]) if pd.notna(sig.get("hs")) else None,
                    "coupling_hz": coupling,
                    "peak_ppm": ppm_list,
                    "peak_heights": heights_list,
                }
            )
        )
    return result


@app.post(
    "/nmr/query",
    response_model=list[NmrQueryResult],
    summary="Query metabolites by 1H-NMR signal pattern",
    tags=["NMR"],
)
def query_nmr(body: NmrQueryRequest):
    """
    Find metabolites whose recorded NMR signals match all provided query
    signals (AND logic).

    Each signal specifies a **chemical shift range** and a **multiplicity**.
    Optionally, peak positions (`ppm`) and relative intensities (`heights`)
    enable multiplet-shape scoring for more precise matching.

    Results are ranked by a similarity score that combines signal coverage
    and multiplet alignment.
    """
    query = [s.model_dump() for s in body.signals]
    # Convert ppm/heights lists to numpy arrays for hmdb_local_tools
    for sig in query:
        if sig.get("ppm") is not None:
            sig["ppm"] = np.array(sig["ppm"], dtype=float)
            sig["heights"] = np.array(sig["heights"], dtype=float)
        else:
            sig.pop("ppm", None)
            sig.pop("heights", None)
        sig["range"] = tuple(sig["range"])

    result_df = multiple_query(query, state.nmrdb, state.df)
    top = result_df.head(body.limit)
    return [
        NmrQueryResult(accession=row["accession"], name=row["name"], similarity=row["similarity"])
        for _, row in top.iterrows()
    ]


# ---------------------------------------------------------------------------
# Routes – Utilities
# ---------------------------------------------------------------------------

@app.get(
    "/formula/parse",
    response_model=AtomCounts,
    summary="Parse atom counts from any chemical formula string",
    tags=["Utilities"],
)
def parse_formula(
    formula: Annotated[str, Query(description="Chemical formula, e.g. C6H12N2O")],
):
    """
    Utility endpoint: parse a raw chemical formula string and return the
    per-element atom counts without looking up the database.
    """
    atoms = count_atoms(formula)
    if not atoms:
        raise HTTPException(status_code=422, detail="Could not parse the provided formula.")
    return AtomCounts(formula=formula, atoms=atoms)


@app.get(
    "/health",
    summary="Health check",
    tags=["Utilities"],
)
def health():
    """Returns database row counts to confirm data is loaded."""
    return {
        "status": "ok",
        "metabolites_loaded": len(state.df),
        "nmr_signals_loaded": len(state.nmrdb),
    }


# ---------------------------------------------------------------------------
# Dev entry-point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)