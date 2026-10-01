import cProfile
import pstats

from experiment_11_eventual_integer_scan import (
    generate_primes,
    generate_gaps,
    scan_fixed_P,
)


def profile_workload():

    primes = generate_primes(2_000_000)
    gaps = generate_gaps(primes)

    scan_fixed_P(
        gaps,
        P=1,
        N_start=1,
        N_end=100_000,
        K=60,
    )


if __name__ == "__main__":

    profiler = cProfile.Profile()

    profiler.enable()

    profile_workload()

    profiler.disable()

    stats = pstats.Stats(profiler)

    stats.sort_stats("cumulative")

    stats.print_stats(30)