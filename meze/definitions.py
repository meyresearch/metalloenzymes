import os
import BioSimSpace as bss
import scipy


ROOT_DIRECTORY = os.path.realpath(os.path.join(os.path.dirname(__file__), ".."))

ADJUST_OPTIONS = ["\tdel: delete perturbations by index",
                  "\tadd: add perturbations in format 'add ligand_1 ligand_2 score'",
                  "\tedit: edit existing LOMAP scores in format 'edit index score'",
                  "\ts: save and continue preparing",
                  "\tq: quit"]

FEMTOSECOND = bss.Units.Time.femtosecond
PICOSECOND = bss.Units.Time.picosecond
NANOSECOND = bss.Units.Time.nanosecond
ANGSTROM = bss.Units.Length.angstrom
KELVIN = bss.Units.Temperature.kelvin
ATM = bss.Units.Pressure.atm

BOLTZMANN_CONSTANT = scipy.constants.Boltzmann
AVOGADROS_NUMBER = scipy.constants.Avogadro

COLOURS = {"PINK": "#D0006F",
           "BLUE": "#0099AB"}