#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30/01/2025

@author: Rolando Armas
Yachay Tech University
Phage Therapy Group
"""

''' Call to protein evolutionary algorithm '''

import sys
import timeit
import multiprocessing as mp
import prot_interface.prot_parserI as parser
from prot_interface.logging_config import setup_logging
import logging
import os

setup_logging()
logger = logging.getLogger(__name__)

def _ensure_spawn_start_method():
    """Ensure multiprocessing uses spawn in an idempotent way."""
    try:
        current = mp.get_start_method(allow_none=True)
        if current != "spawn":
            mp.set_start_method("spawn", force=True)
            logger.info("Multiprocessing start method set to spawn")
    except RuntimeError:
        # Already configured by the active process runtime.
        pass


def _parse_flag(raw_value):
    if isinstance(raw_value, bool):
        return raw_value
    text = str(raw_value).strip()
    if "=" in text:
        text = text.split("=", 1)[1]
    return text.lower() == "true"


def _parse_freq(raw_value):
    text = str(raw_value).strip()
    if "=" in text:
        text = text.split("=", 1)[1]
    return int(text)


def run_eaprot(
    scenario,
    algo_name,
    sim_param_str,
    algo_param_str,
    fitness_idxs_str,
    checkpoint,
    freq,
    mobj,
    randomseed,
    output_root="/output",
):
    _ensure_spawn_start_method()

    tic = timeit.default_timer()
    randomseed = int(randomseed)
    checkpoint = _parse_flag(checkpoint)
    mobj = _parse_flag(mobj)
    freq = _parse_freq(freq)

    algo_params = parser.parse_params(algo_param_str)
    sim_params = parser.parse_params(sim_param_str)
    fitness_idxs = parser.parse_list(fitness_idxs_str)

    output_dir = os.path.join(output_root, "run" + str(randomseed))
    os.makedirs(output_dir, exist_ok=True)
    checkpoint_path = os.path.join(output_dir, "checkpoint.pkl")
    if not os.path.exists(checkpoint_path):
        with open(checkpoint_path, "wb") as f:
            f.write(b"")

    logger.info(
        "Starting eaprot run: scenario=%s algo=%s seed=%s mobj=%s checkpoint=%s",
        scenario, algo_name, randomseed, mobj, checkpoint
    )
    logger.info("Output root: %s | Output dir: %s", output_root, output_dir)
    logger.info("ALGO_PARAMS: %s", algo_params)
    logger.info("SIM_PARAMS: %s", sim_params)
    logger.info("FITNESS_IDXS: %s", fitness_idxs)

    if mobj:
        algo_name_normalized = str(algo_name).strip().lower()

        if algo_name_normalized in {"moea", "moead"}:
            import prot_mob_pymoo as sga_mob
            logger.info("Multiobjective backend selected: MOEA/D (prot_mob_pymoo)")
        elif algo_name_normalized in {"nsga3", "nsgaiii"}:
            import prot_mob_nsga3_pymoo as sga_mob
            logger.info("Multiobjective backend selected: NSGA-III (prot_mob_nsga3_pymoo)")
        else:
            raise ValueError(
                f"Unsupported multiobjective algorithm '{algo_name}'. "
                "Valid options are: moea, moead, nsga3, nsgaiii."
            )

        sga_mob.pymoo_sga_protein(
            scenario, algo_params, sim_params, fitness_idxs, output_dir, randomseed
        ).run(checkpoint=checkpoint, freq=freq)
    else:
        logger.info("mobj=False: no evolutionary execution was launched.")

    toc = timeit.default_timer()
    logger.info("Finished eaprot run: seed=%s elapsed=%.2f seconds", randomseed, toc - tic)
    return output_dir


def main():
    if len(sys.argv) != 10:
        print(len(sys.argv))
        logging.critical(
            "usage: %s <scenario> <sea/moea/nsga3> <sim params> <algo params> "
            "<fitness_idxs> <checkpoint> <freq> <multiobj> <RANDOMSEED>",
            sys.argv[0],
        )
        sys.exit(-1)

    run_eaprot(
        scenario=sys.argv[1],
        algo_name=sys.argv[2],
        sim_param_str=sys.argv[3],
        algo_param_str=sys.argv[4],
        fitness_idxs_str=sys.argv[5],
        checkpoint=sys.argv[6],
        freq=sys.argv[7],
        mobj=sys.argv[8],
        randomseed=sys.argv[9],
    )


if __name__ == "__main__":
    main()
