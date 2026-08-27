import sys
import numba

import pandas as pd
import numpy as np
from numba import jit
import time

# CONFIGURATION

EXCEL_FILE = "patient28_celltypes.xlsx"
WINDOW_SIZE = 700

# Match cell type proportions in the initial spatial transcriptomics data
TARGET_CELL_PROPORTIONS = {
    "CAF": 0.33880,
    "Tumour epithelial": 0.61355,
    "CD8 T cell": 0.04765
}

# Match CD274 expression proportions in the initial spatial transcriptomics data
TARGET_CD274_PROPORTIONS = {
    "Tumour epithelial": 0.071191,  # % of tumour cells expressing CD274
    "CAF": 0.11358,                 # % of CAFs expressing CD274
    "CD8 T cell": 0.088301          # % of CD8 T cells expressing CD274
}

# Weights for balancing the two objectives
CELL_TYPE_WEIGHT = 1.0
CD274_WEIGHT = 1.0

# NUMBA JIT

@jit(nopython=True)
def count_cells_in_window_jit(cells, cd274, x_start, y_start, window_size):
    """
    Count cells by type and CD274+ status in a window
    Returns: [count_CAF, count_Tumour, count_CD8,
              cd274_CAF, cd274_Tumour, cd274_CD8]
    """
    counts = np.zeros(6, dtype=np.int32)
    x_end = x_start + window_size
    y_end = y_start + window_size

    for i in range(cells.shape[0]):
        x = cells[i, 0]
        y = cells[i, 1]
        cell_type = cells[i, 2]
        cd274_val = cd274[i]

        if x_start <= x < x_end and y_start <= y < y_end:
            # Count total cells by type
            counts[cell_type] += 1

            # Count CD274+ cells by type
            if cd274_val > 0:
                counts[3 + cell_type] += 1

    return counts


@jit(nopython=True)
def calculate_score_jit(counts, cell_targets, cd274_targets,
                        cell_type_weight, cd274_weight):
    """
    Calculate score balancing cell type proportions and CD274 expression
    counts: [count_CAF, count_Tumour, count_CD8, cd274_CAF, cd274_Tumour, cd274_CD8]
    Lower score is better
    """
    count_CAF = counts[0]
    count_Tumour = counts[1]
    count_CD8 = counts[2]
    cd274_CAF = counts[3]
    cd274_Tumour = counts[4]
    cd274_CD8 = counts[5]

    total_cells = count_CAF + count_Tumour + count_CD8

    if total_cells == 0:
        return 1e10

    cell_type_score = 0.0 # Cell type proportion error

    actual_caf_prop = count_CAF / total_cells
    actual_tumour_prop = count_Tumour / total_cells
    actual_cd8_prop = count_CD8 / total_cells

    cell_type_score += (actual_caf_prop - cell_targets[0]) ** 2
    cell_type_score += (actual_tumour_prop - cell_targets[1]) ** 2
    cell_type_score += (actual_cd8_prop - cell_targets[2]) ** 2

    cd274_score = 0.0     # CD274 expression error within each cell type

    if count_CAF > 0:
        actual_cd274_caf = cd274_CAF / count_CAF
        cd274_score += (actual_cd274_caf - cd274_targets[0]) ** 2

    if count_Tumour > 0:
        actual_cd274_tumour = cd274_Tumour / count_Tumour
        cd274_score += (actual_cd274_tumour - cd274_targets[1]) ** 2

    if count_CD8 > 0:
        actual_cd274_cd8 = cd274_CD8 / count_CD8
        cd274_score += (actual_cd274_cd8 - cd274_targets[2]) ** 2

    return cell_type_weight * cell_type_score + cd274_weight * cd274_score


# MAIN SCRIPT

def load_and_normalize_data(excel_file):
    """Load Excel file and normalize cell types"""
    print(f"Loading {excel_file}...")
    df = pd.read_excel(excel_file)

    cell_type_col = df.columns[0]
    x_col = df.columns[1]
    y_col = df.columns[2]

    # Normalize Tumour epithelial variants to one type
    df[cell_type_col] = df[cell_type_col].str.replace(
        r"Tumour epithelial \(proliferative\)",
        "Tumour epithelial",
        regex=True
    )

    print(f"Found {len(df)} cells")
    print(f"Cell types: {df[cell_type_col].unique()}")

    # Round coordinates to nearest integer site
    df[x_col] = np.round(df[x_col]).astype(int)
    df[y_col] = np.round(df[y_col]).astype(int)

    x_min, x_max = df[x_col].min(), df[x_col].max()
    y_min, y_max = df[y_col].min(), df[y_col].max()

    print(f"Grid extent: x=[{x_min}, {x_max}], y=[{y_min}, {y_max}]")
    print(f"Grid size: {x_max - x_min + 1} × {y_max - y_min + 1}")

    total_windows = (x_max - x_min - WINDOW_SIZE + 1) * (y_max - y_min - WINDOW_SIZE + 1)
    print(f"Total windows to search: {total_windows:,}")

    # Encode cell types as integers
    cell_type_map = {
        "CAF": 0,
        "Tumour epithelial": 1,
        "CD8 T cell": 2
    }

    df['cell_type_id'] = df[cell_type_col].map(cell_type_map)

    return df, x_col, y_col, cell_type_map, x_min, x_max, y_min, y_max


def calculate_overlap(x1_start, y1_start, x2_start, y2_start, window_size):
    """Calculate overlap percentage between two windows"""
    x1_end = x1_start + window_size
    y1_end = y1_start + window_size
    x2_end = x2_start + window_size
    y2_end = y2_start + window_size

    overlap_x = max(0, min(x1_end, x2_end) - max(x1_start, x2_start))
    overlap_y = max(0, min(y1_end, y2_end) - max(y1_start, y2_start))

    overlap_area = overlap_x * overlap_y
    window_area = window_size * window_size
    overlap_percent = overlap_area / window_area

    return overlap_percent


def find_optimal_windows_exhaustive(cells_array, cd274_array, x_min, x_max, y_min, y_max,
                                     window_size, cell_targets, cd274_targets,
                                     cell_type_weight, cd274_weight, top_n=5, max_overlap=0.5):
    """Find top N non-overlapping windows by exhaustive search"""

    print(f"\nSearching all windows exhaustively...")
    print(f"(Using Numba JIT compilation)")
    print(f"(Selecting top {top_n} windows with max {max_overlap*100:.0f}% overlap)")

    # First pass: find all good windows
    all_windows = []

    start_time = time.time()
    windows_checked = 0

    # Prepare targets as numpy arrays
    cell_targets_arr = np.array(list(cell_targets.values()), dtype=np.float64)
    cd274_targets_arr = np.array(list(cd274_targets.values()), dtype=np.float64)

    total_windows = (x_max - x_min - window_size + 1) * (y_max - y_min - window_size + 1)

    for x_start in range(x_min, x_max - window_size + 1):
        for y_start in range(y_min, y_max - window_size + 1):
            # Count cells using JIT-compiled function
            counts = count_cells_in_window_jit(cells_array, cd274_array, x_start, y_start, window_size)

            # Calculate score using JIT-compiled function
            score = calculate_score_jit(counts, cell_targets_arr, cd274_targets_arr,
                                       cell_type_weight, cd274_weight)

            all_windows.append((score, x_start, y_start, counts.copy()))

            windows_checked += 1

            # Progress indicator
            if windows_checked % 100000 == 0:
                elapsed = time.time() - start_time
                rate = windows_checked / elapsed
                remaining = (total_windows - windows_checked) / rate if rate > 0 else 0
                print(f"  {windows_checked:,}/{total_windows:,} ({rate:.0f} windows/sec, ETA: {remaining/60:.1f} min)")

    elapsed = time.time() - start_time
    print(f"Completed in {elapsed:.1f} seconds ({total_windows/elapsed:.0f} windows/sec)")

    # Second pass: select top N non-overlapping windows
    print(f"\nSelecting non-overlapping windows...")
    all_windows.sort(key=lambda x: x[0])

    selected_windows = []

    for score, x_start, y_start, counts in all_windows:
        # Check if this window overlaps >50% with any already selected window
        overlaps_too_much = False
        for _, sel_x, sel_y, _ in selected_windows:
            overlap = calculate_overlap(x_start, y_start, sel_x, sel_y, window_size)
            if overlap > max_overlap:
                overlaps_too_much = True
                break

        if not overlaps_too_much:
            selected_windows.append((score, x_start, y_start, counts))
            if len(selected_windows) >= top_n:
                break

    return selected_windows


def report_results(top_windows, cell_type_map, cell_targets, cd274_targets, window_size):
    """Print top 5 results as a table"""

    print("\n" + "="*120)
    print(f"TOP 5 NON-OVERLAPPING WINDOWS")
    print("="*120)

    for rank, (score, x_start, y_start, counts) in enumerate(top_windows, 1):
        if counts is None:
            continue

        print(f"\nRank {rank}: x=[{x_start}, {x_start + window_size}), y=[{y_start}, {y_start + window_size})")
        print(f"Score: {score:.6f}")

        # Extract counts
        count_CAF = counts[0]
        count_Tumour = counts[1]
        count_CD8 = counts[2]
        cd274_CAF = counts[3]
        cd274_Tumour = counts[4]
        cd274_CD8 = counts[5]

        total_cells = count_CAF + count_Tumour + count_CD8
        print(f"Total cells: {total_cells}")

        print("\nCell type proportions:")
        print(f"{'Type':<30} {'Target':>12} {'Actual':>12} {'Error':>12}")
        print("-" * 66)

        actual_caf = count_CAF / total_cells if total_cells > 0 else 0
        actual_tumour = count_Tumour / total_cells if total_cells > 0 else 0
        actual_cd8 = count_CD8 / total_cells if total_cells > 0 else 0

        print(f"{'CAF':<30} {cell_targets['CAF']:>11.1%} {actual_caf:>12.1%} {actual_caf - cell_targets['CAF']:>+11.1%}")
        print(f"{'Tumour epithelial':<30} {cell_targets['Tumour epithelial']:>11.1%} {actual_tumour:>12.1%} {actual_tumour - cell_targets['Tumour epithelial']:>+11.1%}")
        print(f"{'CD8 T cell':<30} {cell_targets['CD8 T cell']:>11.1%} {actual_cd8:>12.1%} {actual_cd8 - cell_targets['CD8 T cell']:>+11.1%}")

        print("\nCD274+ expression (within each cell type):")
        print(f"{'Type':<30} {'Target':>12} {'Actual':>12} {'Count':>12}")
        print("-" * 66)

        actual_cd274_caf = cd274_CAF / count_CAF if count_CAF > 0 else 0
        actual_cd274_tumour = cd274_Tumour / count_Tumour if count_Tumour > 0 else 0
        actual_cd274_cd8 = cd274_CD8 / count_CD8 if count_CD8 > 0 else 0

        print(f"{'CAF':<30} {cd274_targets['CAF']:>11.1%} {actual_cd274_caf:>12.1%} {cd274_CAF}/{count_CAF}")
        print(f"{'Tumour epithelial':<30} {cd274_targets['Tumour epithelial']:>11.1%} {actual_cd274_tumour:>12.1%} {cd274_Tumour}/{count_Tumour}")
        print(f"{'CD8 T cell':<30} {cd274_targets['CD8 T cell']:>11.1%} {actual_cd274_cd8:>12.1%} {cd274_CD8}/{count_CD8}")


# RUN

if __name__ == "__main__":
    print("="*120)
    print("Optimal Window Finder (Exhaustive Search with CD274 Expression)")
    print("="*120)

    # Load data
    df, x_col, y_col, cell_type_map, x_min, x_max, y_min, y_max = load_and_normalize_data(EXCEL_FILE)

    # Convert to numpy arrays
    cells_array = df[[x_col, y_col, 'cell_type_id']].values.astype(np.int32)
    cd274_array = df['CD274'].values.astype(np.float32)

    print(f"\nCompiling Numba JIT functions (first run takes ~5-10 seconds)...")

    # Find top 5 non-overlapping windows
    top_windows = find_optimal_windows_exhaustive(
        cells_array, cd274_array, x_min, x_max, y_min, y_max,
        WINDOW_SIZE, TARGET_CELL_PROPORTIONS, TARGET_CD274_PROPORTIONS,
        CELL_TYPE_WEIGHT, CD274_WEIGHT, top_n=5, max_overlap=0.5
    )

    # Report results
    report_results(top_windows, cell_type_map, TARGET_CELL_PROPORTIONS, TARGET_CD274_PROPORTIONS, WINDOW_SIZE)
