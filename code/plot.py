"""Render nine simulation figures and export their plotting data."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results'

def load():
    data = pd.read_csv(DATA / 'summary.csv')
    assert len(data) == 700
    return data.rename(columns={'rate':'rejection_rate'})

COLORS = {'lr': '#0072B2', 'ols': '#D55E00', 'dof': '#009E73', 'loo': '#CC79A7', 'cloo': '#E69F00', 'ploo': '#7B3294', 'oracle': '#333333'}
MARKERS = {'lr': 'o', 'ols': 'o', 'dof': 's', 'loo': 'v', 'cloo': 's', 'ploo': 'D', 'oracle': '^'}
ALL = ['lr', 'ols', 'dof', 'loo', 'cloo', 'ploo', 'oracle']
MAIN = [m for m in ALL if m != 'ols']

def render():
    s = load()
    REC.mkdir(exist_ok=True)
    jobs = [('introduction_CRE', 'CRE', 'normal', ALL, True), ('CRE_main', 'CRE', 'normal', MAIN, False), ('CAR_main', 'CAR', 'normal', MAIN, False), ('BRE_appendix', 'BRE', 'normal', MAIN, False), ('CRE_with_ols', 'CRE', 'normal', ALL, False), ('BRE_with_ols', 'BRE', 'normal', ALL, False), ('CAR_with_ols', 'CAR', 'normal', ALL, False), ('sample_size_appendix', 'CRE', 'normal', ALL, True), ('gamma_CRE_appendix', 'CRE', 'gamma', ['lr', 'dof', 'ploo', 'oracle'], False)]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    points = []
    layouts = []
    for name, design, experiment, methods, one in jobs:
        d = s.loc[s.design.eq(design) & s.experiment.eq(experiment)].copy()
        if name == 'sample_size_appendix':
            d = d.loc[d.ratio.eq('1to1') & d.p.eq(20)]
            xcol = 'n'
        else:
            d = d.loc[d.n.eq(150 if design == 'CAR' else 100)]
            xcol = 'p'
        if name == 'introduction_CRE':
            d = d.loc[d.ratio.eq('1to1')]
        fig, axes = plt.subplots(1 if one else 2, 2, figsize=(10.5, 4.2) if one else (10, 7.4), squeeze=False)
        if one:
            fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.25, wspace=0.25)
        else:
            fig.subplots_adjust(left=0.085, right=0.98, top=0.94, bottom=0.15, hspace=0.38, wspace=0.24)
        for i, ratio in enumerate(['1to1'] if one else ['1to1', '2to1']):
            for j, scope in enumerate(['H0', 'H1']):
                ax = axes[i, j]
                levels = sorted(d[xcol].unique())
                if j == 0:
                    band = 0.025 + np.array([-1, 1]) * 1.95996398454 * np.sqrt(0.025 * 0.975 / 10000)
                    ax.axhspan(*band, color='#F0F0F0')
                    ax.axhline(0.025, color='#777777', linewidth=0.9, linestyle='-')
                for method in methods:
                    dd = d.loc[d.ratio.eq(ratio) & d.scope.eq(scope) & d.method.eq(method)].sort_values(xcol)
                    assert dd[xcol].tolist() == levels and (not dd[xcol].duplicated().any())
                    label = f'$\\hat{{\\chi}}_{{\\mathrm{{{method}}}}}$'
                    ax.plot(np.arange(len(levels)), dd.rejection_rate, color=COLORS[method], marker=MARKERS[method], linestyle='-', linewidth=1.8, markersize=4.5, label=label)
                    for r in dd.itertuples():
                        points.append(dict(figure=name, panel=f'{ratio}_{scope}', group=r.group, p=int(r.p), n=int(r.n), method=method, value=r.rejection_rate, valid=int(r.valid), invalid=int(r.invalid)))
                ax.set_xlabel('$p$' if xcol == 'p' else '$n$')
                ax.set_ylabel('Type I error' if j == 0 else 'Power')
                ax.set_xticks(np.arange(len(levels)), [str(int(v)) for v in levels])
                ax.set_xlim(-0.18, len(levels) - 0.82)
                top = max(d.loc[d.scope.eq(scope) & d.method.isin(methods), 'rejection_rate'].max() * 1.1, 0.04 if j == 0 else 0.1)
                ax.set_ylim(0, top)
                ax.grid(axis='y', color='#E3E3E3', linewidth=0.6)
                ax.spines[['top', 'right']].set_visible(False)
                ax.set_title(f'$H_{{{j}}}$' if one else f"{('1:1' if i == 0 else '2:1')}: $H_{{{j}}}$")
        handles, labels = axes.flat[0].get_legend_handles_labels()
        leg = fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.52, 0.025), ncol=len(methods), frameon=False, handlelength=1.7, columnspacing=1.1)
        fig.canvas.draw()
        ren = fig.canvas.get_renderer()
        boxes = [ax.get_tightbbox(ren) for ax in axes.flat] + [leg.get_window_extent(ren)]
        for i, b in enumerate(boxes):
            assert b.x0 >= 0 and b.y0 >= 0 and (b.x1 <= fig.bbox.width) and (b.y1 <= fig.bbox.height), (name, 'bounds')
            assert all((not b.overlaps(c) for c in boxes[i + 1:])), (name, 'overlap')
        fig.savefig(FIG / f'{name}.pdf')
        fig.savefig(FIG / f'{name}.png', dpi=180)
        plt.close(fig)
        layouts.append(dict(figure=name, methods=methods, equal_spacing=True, solid_lines=True, layout_passed=True))
    p = pd.DataFrame(points)
    v = p.merge(s[['group', 'p', 'method', 'rejection_rate', 'valid', 'invalid']], on=['group', 'p', 'method'], validate='many_to_one', suffixes=('', '_raw'))
    assert len(p) == 1134 and np.max(abs(v.value - v.rejection_rate)) < 1e-15
    assert v.valid.eq(v.valid_raw).all() and v.invalid.eq(v.invalid_raw).all()
    
    for name, values in p.groupby('figure', sort=False):
        values.to_csv(REC / (name + '.csv'), index=False)
    (FIG.parent / 'plot_check.json').write_text(json.dumps(layouts, indent=2))
    print('Saved nine figures and their plotting data.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/figures')
    args = parser.parse_args()
    FIG = args.output.resolve()
    REC = FIG.parent / 'plot_data'
    FIG.mkdir(parents=True, exist_ok=True)
    REC.mkdir(parents=True, exist_ok=True)
    render()
