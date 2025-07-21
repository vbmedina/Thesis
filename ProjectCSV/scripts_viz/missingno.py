def compute_matrix_completeness(df):
    """
    Compute matrix completeness percentage. Either provide a DataFrame, or provide the components directly:
    - unique_cells: list/set or int
    - unique_compounds: list/set or int
    - data_points: number of measurements
    :return: Completeness percentage (0-100)
    """
    #number of na cells
    total_missing = df.isna().values.sum()
    return 100 - (total_missing / (df.shape[0] * df.shape[1]) * 100)

