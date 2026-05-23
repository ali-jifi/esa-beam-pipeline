PROBES = ("a", "b", "c", "d", "e")


def eflux_var(probe: str, species: str = "i") -> str:
    """tha_peif_en_eflux for probe='a', species='i'."""
    return f"th{probe}_pe{species}f_en_eflux"
