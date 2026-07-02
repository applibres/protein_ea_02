#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sequential replicate runner for eaprot.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time

from eaprot_call import run_eaprot
from prot_interface.logging_config import setup_logging


setup_logging()
logger = logging.getLogger(__name__)


def _parse_bool(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {value}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run a sequence of eaprot replicas and isolate outputs by replicate set."
    )
    parser.add_argument("--init", type=int, required=True, help="Initial random seed (inclusive).")
    parser.add_argument("--end", type=int, required=True, help="Final random seed (inclusive).")
    parser.add_argument("--scenario", required=True, help="Scenario name.")
    parser.add_argument("--algo", required=True, help="Algorithm name: sea/moea.")
    parser.add_argument("--sim-params", required=True, help="Simulation params string.")
    parser.add_argument("--algo-params", required=True, help="Algorithm params string.")
    parser.add_argument("--fitness-idxs", required=True, help="Fitness list string, e.g. fitness_idxs=0,2,7.")
    parser.add_argument("--checkpoint", type=_parse_bool, required=True, help="Restore from checkpoint.")
    parser.add_argument("--freq", type=int, required=True, help="Checkpoint frequency.")
    parser.add_argument("--mobj", type=_parse_bool, required=True, help="Use multiobjective mode.")
    parser.add_argument(
        "--replicates-id",
        type=int,
        default=1,
        help="Replica-set namespace id. Output is stored in /output/replicates#.",
    )
    parser.add_argument(
        "--sleep-seconds",
        type=int,
        default=30,
        help="Sleep time between runs (seconds).",
    )
    parser.add_argument(
        "--output-root",
        default="/output",
        help="Base output root (default: /output).",
    )
    return parser


def main():
    args = build_parser().parse_args()

    if args.replicates_id < 1:
        logger.error("replicates-id must be >= 1")
        return 1
    if args.end < args.init:
        logger.error("end must be >= init")
        return 1
    if args.sleep_seconds < 0:
        logger.error("sleep-seconds must be >= 0")
        return 1

    tic = time.time()
    output_root = os.path.join(args.output_root, f"replicates{args.replicates_id}")
    logger.info(
        "Starting replicate block: seeds=[%s,%s] scenario=%s mobj=%s output_root=%s",
        args.init, args.end, args.scenario, args.mobj, output_root
    )

    for seed in range(args.init, args.end + 1):
        run_tic = time.time()
        logger.info("Starting replicate seed=%s", seed)
        try:
            run_dir = run_eaprot(
                scenario=args.scenario,
                algo_name=args.algo,
                sim_param_str=args.sim_params,
                algo_param_str=args.algo_params,
                fitness_idxs_str=args.fitness_idxs,
                checkpoint=args.checkpoint,
                freq=args.freq,
                mobj=args.mobj,
                randomseed=seed,
                output_root=output_root,
            )
        except Exception:
            logger.exception("Error in replicate seed=%s. Aborting replicate block.", seed)
            return 1

        logger.info(
            "Finished replicate seed=%s output=%s elapsed=%.2f seconds",
            seed, run_dir, time.time() - run_tic
        )

        if seed < args.end and args.sleep_seconds > 0:
            time.sleep(args.sleep_seconds)

    total = time.time() - tic
    logger.info("Successful replicate block. Total elapsed: %.2f seconds", total)
    return 0


if __name__ == "__main__":
    sys.exit(main())


'''
python replicates.py \
  --init 1 --end 1 \
  --scenario test06 \
  --algo moea \
  --sim-params "pdbfile=9Q1V_prepared_clean_relaxed.pdb,partners=A_B,ligand_chain=B" \
  --algo-params "gen=20,popsize=25,mutp=0.3,bo_enabled=1,bo_candidates_per_parent=8,bo_beta=1.0,bo_min_train=24" \
  --fitness-idxs "fitness_idxs=7" \
  --checkpoint false \
  --freq 2 \
  --mobj true \
  --replicates-id 1
'''
