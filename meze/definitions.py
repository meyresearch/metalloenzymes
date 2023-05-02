import os

ROOT_DIRECTORY = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))

ADJUST_OPTIONS = ["\tdel: delete perturbations by index",
                  "\tadd: add perturbations in format 'add ligand_1 ligand_2 score'",
                  "\tedit: edit existing LOMAP scores in format 'edit index score'",
                  "\ts: save and continue preparing",
                  "\tq: quit"]