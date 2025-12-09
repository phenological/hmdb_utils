from rapidfuzz import process, fuzz
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import ast
import re

import bisect

def approximate_lookup(df, column, search_string,
                       scorer1=fuzz.partial_token_ratio,
                       scorer2=fuzz.ratio,
                       limit=10):
    """
    Perform approximate string matching lookup on a specified column of a DataFrame.

    This function searches for entries similar to a given search string within
    a specified column. It applies two sequential scoring functions:

    1. **scorer1** (default: `fuzz.partial_token_ratio`)  
       Used to retrieve all candidate rows from the target column using
       `rapidfuzz.process.extract`.

    2. **scorer2** (default: `fuzz.ratio`)  
       Applied to each candidate’s main name and its synonyms to determine the
       best-matching string.

    The function returns the top results ranked by the secondary score.

    Args:
        df (pd.DataFrame):
            DataFrame containing at least the columns:
            - the target search column,
            - `"name"` (main identifier),
            - `"synonyms"` (iterable of synonyms for the row).
        column (str):
            Name of the column in `df` to perform the lookup on.
        search_string (str):
            The string to search for.
        scorer1 (Callable, optional):
            A `rapidfuzz` scoring function used for initial candidate extraction.
            Defaults to `fuzz.partial_token_ratio`.
        scorer2 (Callable, optional):
            A `rapidfuzz` scoring function used to score the main name and its
            synonyms for each candidate. Defaults to `fuzz.ratio`.
        limit (int, optional):
            Maximum number of results to return. Defaults to 10.

    Returns:
        np.ndarray or None:
            A NumPy array with shape `(n, 4)` containing:
            `[main_name, best_synonym, row_index, score]`,  
            sorted by score in descending order.

            Returns `None` if no matches are found.

    Notes:
        - Values in the search column and the search string are lowercased
          before matching.
        - The `"synonyms"` column must contain an iterable (e.g., list) of
          synonym strings. If it contains a single comma-separated string,
          it should be split beforehand.


    Example:
        >>> approximate_lookup(df, 'synonyms_cat', 'citric acid')
    """

    # Extract the column values as a list
    choices = df[column].str.lower().tolist()
    names = df['name'].values
    search_string = search_string.lower()
    # Get the best matches
    matches = process.extract(search_string, choices, scorer=scorer1, score_cutoff=100, limit=None)
    #print(matches)
    if len(matches) > 0:
        result = []
        for match in matches:
            matchings_names = df.loc[match[2],['synonyms']].values[0]
            best_name = ""
            best_score = 0
            for name in [names[match[2]]] + matchings_names:
                score = scorer2(search_string, name.lower())
                if score > best_score:
                    best_score = score
                    best_name = name
            result.append([names[match[2]], best_name, match[2], best_score,])
        # Map to the main name
        matches = pd.DataFrame(result).sort_values(3, ascending=False)
        #for synomym in data['synomyns'].
        return matches.head(limit).values
    return None


def triangle_overlap_area(x_a, y_a, x_b, y_b, b):
    """
    Compute the overlap area between two isosceles triangles with:
        - same base length b
        - bases on y = 0
        - centers at x_a and x_b on the x-axis
        - heights y_a and y_b
    """
    def tri_height(x, xc, h, b):
        """Return the height of the triangle at position x (piecewise linear)."""
        half = b / 2
        left, right = xc - half, xc + half
        if x < left or x > right:
            return 0.0
        # Linear interpolation
        if x <= xc:
            return h * (x - left) / (xc - left)
        else:
            return h * (right - x) / (right - xc)

    # Compute global horizontal overlap range
    L = max(x_a - b/2, x_b - b/2)
    R = min(x_a + b/2, x_b + b/2)

    if R <= L:
        return 0.0  # No overlap

    # Integrate numerically (analytical solution is possible but long)
    # Fine grid for accuracy
    xs = np.linspace(L, R, 2001)
    h_min = np.minimum(
        [tri_height(x, x_a, y_a, b) for x in xs],
        [tri_height(x, x_b, y_b, b) for x in xs]
    )

    # Numerical integration using trapezoid rule
    area = np.trapz(h_min, xs)
    return area


# Range query function
def _range_query_and(query, nmrdb):
    indexer = [(a, b) for a, b in zip(nmrdb['shift1(ppm)'].values, nmrdb.index.values)]
    result_ids = None
    for rng in query:
        range_i = rng['range']
        left_index = bisect.bisect_left(indexer,  (range_i[0], -1))
        right_index = bisect.bisect_right(indexer, (range_i[1], float('inf')))
        res = set(nmrdb.loc[[item[1] for item in indexer[left_index:(right_index + 1)]],['accession']]['accession'].values)
        if result_ids is None:
            result_ids = res
        else:
            result_ids = result_ids & res
    return list(result_ids)

def _multiplet_match(test_ppm, test_heights, query_ppm, query_heights):
    #1.  Center both multiplets at center of mass
    cm_test = np.sum(test_ppm * test_heights) / np.sum(test_heights)
    test_ppm -= cm_test
    cm_query = np.sum(query_ppm * query_heights) / np.sum(query_heights)
    query_ppm -= cm_query
    #2. Compare area overlap when both multiplets has area 1.
    score =0
    base = 0.008 # TODO: We could change this. 
    for i in range(test_ppm.shape[0]):
        for j in range(query_ppm.shape[0]):
            score += triangle_overlap_area(test_ppm[i], test_heights[i], query_ppm[j], query_heights[j], base)
    assert (score >=0) & (score <= 1), "Score out of range"
    return score * 2 / base    
    

def multiple_query(query, nmrdb, metabolytes):
    """
    Perform multi-signal lookup in an HMDB-like 1H NMR peak database.

    This function allows querying metabolites based on one or more expected NMR
    signals. Each signal may specify:

    - a **chemical shift range** where the peak must occur,
    - a **multiplicity** to match (e.g., "d", "t", "q", or "*" for any),
    - optionally a **multiplet shape** defined by peak positions (`ppm`) and
      corresponding intensities (`heights`).

    If multiplet shape information is provided, the function applies an
    additional refinement step by computing multiplet similarity through
    `_multiplet_match`, allowing more accurate matching of complex peak patterns.

    The similarity score returned for each metabolite combines:
    - the proportion of expected signals found, normalized by the number of
      signals present for that metabolite, and
    - bonus similarity for matching multiplet shapes.

    Args:
        query (list of dict):
            A list of signal descriptors. Each dictionary must include:
            - `range` (tuple): `(from_ppm, to_ppm)` search interval.
            - `mult` (str): Peak multiplicity to match, or `"*"` for wildcard.

            Optionally, each dictionary may also include:
            - `ppm` (list or array-like): Ideal multiplet peak positions.
            - `heights` (list or array-like): Relative intensities for the peaks.
              These will be normalized internally.

            Example element:
            `{'range': (4.0, 4.2), 'mult': 'q', 'ppm': [...], 'heights': [...]}`

        nmrdb (pd.DataFrame):
            A long-format peak table containing at least:
            - `"accession"` (metabolite identifier),
            - `"name"` (metabolite name),
            - `"shift1(ppm)"` (peak position),
            - `"type"` (multiplicity),
            - For multiplet matching: `"ppm"` and `"heights"` per signal.

        metabolytes (pd.DataFrame):
            A metabolite annotation table with at least the `"accession"` column.
            Used to filter results after peak-level querying.

    Returns:
        pd.DataFrame:
            A DataFrame sorted by descending similarity with columns:
            - `"accession"`: Metabolite identifier.
            - `"name"`: Metabolite name.
            - `"similarity"`: Computed match score based on:
                * number of matching signals,
                * normalized multiplet similarity when applicable.

    Notes:
        - If a query signal contains `ppm` and `heights`, the algorithm attempts
          to compute multiplet alignment via `_multiplet_match` after confirming
          at least one peak exists in the defined range.
        - Intensities (`heights`) are normalized automatically.
        - The function performs an AND-style query across all provided ranges:
          metabolites must match every query signal at least once to receive a
          high score.
        - The lookup is tolerant: matching multiplicity can be forced with `"*"`.

    Examples:
        # Searching for Leucine
        query = [{'range': (0.94, 0.99), 'mult': 't', 'ppm': np.array([0.949542,0.96010, 0.970836]), 'heights': np.array([0.25,0.5,0.25])},]

        result = multiple_query(query, nmrdb, df)
        result.head(6) 
    """
    res = _range_query_and(query, nmrdb)
    filtered_df = metabolytes.loc[metabolytes['accession'].isin(res), ]
    result = []
    
    # Ensure correct query structure
    for query_i in query:
        if 'ppm' in query_i:
            query_i['ppm'] = np.array(query_i['ppm'])
            query_i['heights'] = np.array(query_i['heights'])
            query_i['heights'] /= np.sum(query_i['heights'])
        
    for row in filtered_df.iterrows():
        signals = nmrdb[nmrdb['accession'] == row[1]["accession"]]
        matches = 0
        for query_i in query:
            range_i = query_i['range']
            mult = query_i['mult']
            from_i = range_i[0]
            to_i = range_i[1]
            mul_i = mult
            for signal in signals.iterrows():
                if signal[1]['shift1(ppm)'] >= from_i and signal[1]['shift1(ppm)'] <= to_i and ((signal[1]['type'] == mul_i) or (mul_i == "*")):
                    matches += 1
                    # Check multiplet alignment
                    if 'ppm' in query_i:
                        matches += _multiplet_match(signal[1]['ppm'], signal[1]['heights'], query_i['ppm'], query_i['heights'])
                    break
        result.append((row[1]["accession"], row[1]["name"], len(query) / signals.shape[0] + matches /  len(query)))
                    
    result = pd.DataFrame(result, columns=["accession", "name", "similarity"])
    result = result.sort_values("similarity", ascending=False)
    return result

def count_atoms(formula):
    """
    Counts the number of atoms of each element in a chemical formula.

    Args:
        formula: A string representing the chemical formula.

    Returns:
        A dictionary where keys are element symbols and values are the number of atoms of that element.
    
    Examples:
        count_atoms("C6H12N2")
    """
    atoms = {}
    for match in re.findall(r"([A-Z][a-z]*)(\d*)", formula):
        element = match[0]
        count = int(match[1]) if match[1] else 1
        atoms[element] = atoms.get(element, 0) + count
    return atoms

# Apply the function to the 'chemical_formula' column
#atom_counts = df['chemical_formula'].apply(count_atoms)

# Expand the dictionary of atom counts into separate columns
#df = pd.concat([df, atom_counts.apply(pd.Series)], axis=1)
#print(df)


    