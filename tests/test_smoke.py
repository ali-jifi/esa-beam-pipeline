from esa_plotting import PROBES, eflux_var


def test_probes_tuple():
    assert PROBES == ("a", "b", "c", "d", "e")


def test_eflux_var_default_species():
    assert eflux_var("a") == "tha_peif_en_eflux"


def test_eflux_var_electron():
    assert eflux_var("c", species="e") == "thc_peef_en_eflux"


def test_package_imports():
    import esa_plotting  # noqa: F401
    from esa_plotting import configure_eflux_panel, load_esa, set_data_dir  # noqa: F401
