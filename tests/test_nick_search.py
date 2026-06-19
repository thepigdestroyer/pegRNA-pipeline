"""Task 5: enumerate candidate nick sites in [40,70] bp on the opposite strand."""

from pegrna_pipeline.hpc.nick_search import NickCandidate, enumerate_nick_sites


def test_enumerate_nick_sites_in_window_opposite_strand():
    # 200-nt genome; pegRNA on Rv => nicks wanted on Fw; S at 100
    genome = "ACGT" * 50
    sites = enumerate_nick_sites(
        genome, site_s=100, pegrna_strand="Rv",
        distance_min=40, distance_max=70,
    )
    assert all(isinstance(s, NickCandidate) for s in sites)
    assert all(40 <= s.distance_to_edit <= 70 for s in sites)
    assert all(s.target_strand == "Fw" for s in sites)  # opposite of Rv pegRNA
    assert all(len(s.protospacer) == 20 for s in sites)
    assert all(s.pam[1:] == "GG" for s in sites)  # NGG


def test_no_pam_in_window_returns_empty():
    genome = "AAAAAAAA" * 30  # no NGG/CCN PAM anywhere
    assert enumerate_nick_sites(
        genome, site_s=120, pegrna_strand="Fw",
        distance_min=40, distance_max=70,
    ) == []


def test_forward_pegrna_wants_reverse_nicks():
    genome = "ACGT" * 50
    sites = enumerate_nick_sites(
        genome, site_s=100, pegrna_strand="Fw",
        distance_min=40, distance_max=70,
    )
    assert all(s.target_strand == "Rv" for s in sites)


def test_distance_window_respected():
    genome = "GGACGTACGTACGT" * 15
    sites = enumerate_nick_sites(
        genome, site_s=100, pegrna_strand="Rv",
        distance_min=45, distance_max=55,
    )
    assert all(45 <= s.distance_to_edit <= 55 for s in sites)
