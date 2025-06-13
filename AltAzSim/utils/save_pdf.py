from AltAzSim.coveragemap.coveragemap import CoverageMap
from AltAzSim.coveragemap.overlay import Overlay

import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
import os


def save_maps_pdf(cv_maps: CoverageMap | list[CoverageMap], save_loc: str, metadata: str | None = None, overlays: Overlay | list[Overlay] | None = None):
    if not isinstance(cv_maps, list):
        cv_maps = [cv_maps]
    if not all(isinstance(cv_map, CoverageMap) for cv_map in cv_maps):
        raise ValueError("all elements of cv_maps must be CoverageMap instances")
    if overlays is not None and not isinstance(overlays, list):
        overlays = [overlays]
    if overlays is not None and not all(isinstance(overlay, Overlay) for overlay in overlays):
        raise ValueError("all elements of overlays must be Overlay instances")
    
    if not save_loc.lower().endswith(".pdf"):
        save_loc += ".pdf"

    # creates folders if needed
    folder = os.path.dirname(save_loc)
    if folder:
        os.makedirs(folder, exist_ok=True)

    with PdfPages(save_loc) as pdf:
        for i, cv_map in enumerate(cv_maps):
            cv_map._prep_for_disp()
            fig = plt.gcf()

            if overlays is not None:
                for overlay in overlays:
                    overlay._place_overlay()

            pdf.savefig(fig, dpi=300, bbox_inches='tight')
            plt.close(fig)

        if metadata is not None:
            _add_text_pages(pdf, metadata)


def _add_text_pages(pdf, text, font_size=8):
    fig_width, fig_height = 8.5, 11  # inches
    line_height = font_size / 72     # points to inches
    max_lines_per_page = int(fig_height / line_height) - 5  # margin buffer

    lines = text.replace('\t', '    ').splitlines()

    for page_start in range(0, len(lines), max_lines_per_page):
        fig = plt.figure(figsize=(fig_width, fig_height))
        plt.axis('off')

        page_lines = lines[page_start:page_start + max_lines_per_page]
        for i, line in enumerate(page_lines):
            y_pos = 1 - (i + 1) * line_height / fig_height
            fig.text(0.05, y_pos, line, fontsize=font_size, family='monospace', ha='left', va='top')

        pdf.savefig(fig)
        plt.close(fig)