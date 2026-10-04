"""Figure saving and frozen table reads; no scientific calculations on import."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from ..run_io import read_table
GROUP_COLOURS = ['#4477aa', '#ee6677', '#228833', '#ccbb44', '#66ccee', '#aa3377', '#bbbbbb', '#332288', '#44aa99']
QC_COLOURS = {'pass': '#228833', 'fail': '#cc3311', 'unknown': '#999999'}

def frame(title, rows=1, columns=1, size=(12, 6)):
    figure = Figure(figsize=size, layout='constrained')
    FigureCanvasAgg(figure)
    axes = figure.subplots(rows, columns, squeeze=False).ravel()
    figure.suptitle(title)
    return (figure, axes)

def save(figure, output, name, caption):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    figure.savefig(output / (name + '.png'), dpi=180)
    figure.savefig(output / (name + '.pdf'))
    (output / (name + '.txt')).write_text(caption + '\n')
    figure.clear()

def aligned_annotations(source):
    ids = (source / 'rmsd_entry_index.txt').read_text().splitlines()
    table = read_table(source, 'entry_annotations.csv')
    if table.chain_key.duplicated().any():
        raise ValueError('Selected annotation keys must be unique')
    return table.set_index('chain_key').loc[ids].reset_index()
