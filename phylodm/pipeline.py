from phylodm.pdm import PDM


def newick_to_pdm(path_newick: str, method: str, path_out: str):
    """Converts a newick file to a PDM and caches it to disk."""
    pdm = PDM.get_from_newick_file(path_newick, method)
    pdm.save_to_path(path_out)
