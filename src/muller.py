import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def _get_strains_ordering(adjacency_df):

    children_by_parent = adjacency_df.groupby('Parent')['Identity'].apply(lambda x: list(sorted(x)))

    def get_inner_order(identity):

        children_identities = children_by_parent.get(identity, [])

        if len(children_identities) == 0:
            return [identity, identity]

        inner = [get_inner_order(c) for c in children_identities]
        return [identity] + sum(inner, []) + [identity]

    order = []

    identities = list(set(adjacency_df['Identity'].values) | set(adjacency_df['Parent'].values))

    for strain in sorted(identities):
        if strain not in order:
            order += get_inner_order(strain)

    return np.array(order)


def _get_y_values(populations_df, adjacency_df, smoothing_std):

    ordering = _get_strains_ordering(adjacency_df)
    population_size_max = populations_df.groupby('Generation')['Population'].sum().max()
    generations = populations_df['Generation'].max() - populations_df['Generation'].min()

    pivot = populations_df.pivot(index='Generation', columns='Identity', values='Population').sort_index()
    pivot = pivot.rolling(generations, 1, True, 'gaussian').mean(std=smoothing_std).clip(0, population_size_max)

    Y = pivot[ordering] / 2

    # Avoid middle lines. Double leaf clones.
    keep = [0]

    for i, c in enumerate(Y.columns[1:], 1):
        if c == Y.columns[i - 1]:
            Y.iloc[:, i] *= 2
            keep.pop()

        keep.append(i)

    return Y.iloc[:, keep]


import matplotlib.pyplot as plt
import matplotlib
import numpy as np

def _muller_plot(populations_df, adjacency_df, color_by, title='',colormap='terrain', colorbar=False, background_strain=True,
                 smoothing_std=10, normalize=False, norm_range=None, ax=None, click_callback=None):

    if normalize:
        populations_df['Population'] = populations_df.groupby('Generation')['Population'].transform(lambda x: x / x.sum())

    x = populations_df['Generation'].unique()
    y_table = _get_y_values(populations_df, adjacency_df, smoothing_std)

    final_order = y_table.columns.values
    Y = y_table.to_numpy().T

    cmap = plt.get_cmap(colormap)
    color_by = color_by.copy()

    final_order_str = [str(i) for i in final_order]
    ordered_colors = color_by.loc[final_order_str]
  
    if normalize and norm_range is None:
            norm_range = (np.min(ordered_colors), np.max(ordered_colors))

    norm = matplotlib.colors.Normalize(vmin=norm_range[0], vmax=norm_range[1])
    colors = cmap(norm(ordered_colors.values))

    if background_strain:
        colors[0] = colors[-1] = [1, 1, 1, 1]

    if ax is None:
        fig, ax = plt.subplots(figsize=(5.46, 5.46*0.25))

    else:
        fig = ax.figure if hasattr(ax, 'figure') else ax.fig

    ax.stackplot(x, Y, colors=colors)

    if colorbar:
        cax = fig.add_axes([0.92, 0.13, 0.02, 0.7])
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
        cb.set_label(color_by.name,rotation=270,labelpad=9)
        cb.ax.yaxis.label.set_fontproperties(matplotlib.font_manager.FontProperties(family='Arial', size=9))
    
    ax.set_xlabel('Time', fontdict={'family': 'Arial', 'size': 9})
    ax.set_ylabel('Frequency' if normalize else 'Abundance', fontdict={'family': 'Arial', 'size': 9})

    plt.sca(ax)

    if click_callback is not None:

        def onclick(event):

            generation = int(np.round(event.xdata))
            row = y_table.loc[generation]
            identity = row[row.cumsum() > event.ydata].index[0]

            click_callback(generation, identity)

        fig.canvas.mpl_connect('button_press_event', onclick)


   # Set additional properties to the exported figure
    ax.tick_params(axis='both', which='major', labelsize=8, labelcolor='black', pad=3)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname('Arial')

    ax.set_yticks([0.0,1.0])
    ax.set_yticks([0.25, 0.5, 0.75], minor=True)
    ax.set_xlabel('Time', fontdict={'family': 'Arial', 'size': 9})
    ax.set_ylabel('Variant frequency', fontdict={'family': 'Arial', 'size': 9})
    ax.set_xlim([0,1001])
    ax.set_ylim([0,1])




    return fig, ax
