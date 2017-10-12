"""
plotting function

:author: Benjamin Schubert
"""
import plotly
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines


def interactive_scatter(df, x_name, y_name, ann_name, color="blue", x_title="", y_title=""):
    """
    generates an interactive scatter plot using plotly in ipython offline mode

    :param df: pandas dataframe which contains all relevant information
    :param x_name: the x column name
    :param y_name: the y column name
    :param ann_name: the annotation column name
    :param color: specifies the color of the dots (default: blue)
    :param x_title: x axis title
    :param y_title: y axis title
    :return: plotly image
    """

    fig = {
        'data': [
            {
                'marker': {'color': color},
                'x': df[x_name],
                'y': df[y_name],
                'text': df[ann_name],
                'mode': 'markers',
            }

        ],
        'layout': {
            'xaxis': {'title': x_title, },
            'yaxis': {'title': y_title}
        }
    }

    # IPython notebook
    return plotly.offline.iplot(fig, filename='pandas/scatter')


def joint_plot(df, x_name, y_name, color="#86a7c5", x_title="", y_title="", out=None):
    """
    generates a joint scatter plot using seaborn

    :param df: pandas dataframe which contains all relevant information
    :param x_name: x column name
    :param y_name: y column name
    :param color: color in hex code
    :param x_title: x axis title
    :param y_title: y axis title
    :param out: output path
    :return: matplotlib figure
    """
    g = sns.jointplot(x_name, y_name, data=df, kind="scatter", color=color, space=0, size=7)
    g.set_axis_labels(x_title, y_title)
    if out is not None:
        g.savefig(out)
    return g


def rank_change_plot(df, label, rank_columns, max_plot=5, rank_names=None, out=None, fontsize=15, cmap=None):
    """
        generates a rank-change plot between specified columns on which ranks are calculated

    :param df: pandas dataframe containing all data
    :param label: the label column
    :param rank_columns: a list of columns which will be plotted
    :param rank_names: optional list of column names which will be plotted instead of the rank_columns
    :param out: optinal path to output file
    :return: matplotlib figure
    """
    def slope_from_points(point1, point2):
        return (point2[1] - point1[1]) / (point2[0] - point1[0])

    def plot_secant(point1, point2, ax):
        # plot the secant
        slope = slope_from_points(point1, point2)
        intercept = point1[1] - slope * point1[0]
        # update the points to be on the axes limits
        x = ax.get_xlim()
        y = ax.get_ylim()
        data_y = [x[0] * slope + intercept, x[1] * slope + intercept]
        line = mlines.Line2D(x, data_y, color='grey')
        ax.add_line(line)

    if cmap is None:
        cmap = sns.color_palette("Set2", max_plot)

    df = df.copy()
    n = len(rank_columns)
    #plt.style.use("seaborn-whitegrid")
    fig = plt.figure(frameon=False)
    plt.box(on=None)
    gs = gridspec.GridSpec(2, n+(n - 1), height_ratios=[1, 7], wspace=0.01, hspace=0.01)

    if rank_names is None or len(rank_names) != n:
        rank_names = rank_columns

    # plot starts here

    # plot first rank
    # first plot distribution and name
    ax11 = plt.subplot(gs[0, 0])
    ax11.hist(df[rank_columns[0]].values, 50, normed=1, alpha=0.70, color=cmap[0])
    ax11.set_title(rank_names[0], fontsize=fontsize)
    ax11.get_yaxis().set_visible(False)
    ax11.tick_params(direction='in')
    ax11.axis('off')
    ax11.set_xlim(0, max(df[rank_columns[0]].values))

    # blot vertical barplot
    ax12 = plt.subplot(gs[1, 0])
    df[rank_columns[0] + "_rank"] = df[rank_columns[0]].rank(ascending=0)
    df.sort_values(by=rank_columns[0] + "_rank", inplace=True)
    label_pref = df.head(max_plot)[label]

    y_pos = [i+0.5 for i in np.arange(len(label_pref))]
    ax12.barh(y_pos, df.head(max_plot)[rank_columns[0]].tolist(), align='center', color=cmap[0])
    ax12.set_yticks(y_pos)
    ax12.set_yticklabels(label_pref)
    ax12.set_ylim(0, max_plot)
    ax12.invert_yaxis()  # labels read top-to-bottom
    ax12.spines['right'].set_visible(False)
    ax12.spines['top'].set_visible(False)
    ax12.grid(False)

    for i in range(1, n):
        j = 2*i # the barplot index
        k = j - 1 # the line plot index

        # plot rank change
        # first plot distribution and name
        ax11 = plt.subplot(gs[0, j])
        ax11.hist(df[rank_columns[i]].values, 50, normed=1, alpha=0.70, color=cmap[i])
        ax11.set_title(rank_names[i], fontsize=fontsize)
        ax11.get_yaxis().set_visible(False)
        ax11.tick_params(direction='in')
        ax11.axis('off')
        ax11.set_xlim(0, max(df[rank_columns[i]].values))

        # blot vertical barplot
        ax12 = plt.subplot(gs[1, j])
        df[rank_columns[i] + "_rank"] = df[rank_columns[i]].rank(ascending=0)
        df.sort_values(by=rank_columns[i] + "_rank", inplace=True)
        ax12.barh(y_pos, df.head(max_plot)[rank_columns[i]].tolist(), align='center', color=cmap[i])
        ax12.set_yticks(y_pos)
        ax12.set_ylim(0, max_plot)
        ax12.invert_yaxis()  # labels read top-to-bottom
        ax12.get_yaxis().set_visible(False)
        ax12.spines['right'].set_visible(False)
        ax12.spines['top'].set_visible(False)
        ax12.grid(False)

        # plot line plot between different ranks
        # inefficient should be done after all ranks are calculated than one one go through is needed.
        ax22 = plt.subplot(gs[1, k])
        for l in label_pref:
            row = df[df[label] == l]

            y1 = row[rank_columns[i - 1] + "_rank"].values - 0.5
            y2 = row[rank_columns[i] + "_rank"].values - 0.5

            plot_secant([0, y1], [1, y2], ax22)

            ax22.set_xlim((0, 1))
            ax22.set_ylim((0, max_plot))
            ax22.invert_yaxis()
            ax22.get_yaxis().set_visible(False)
            ax22.get_xaxis().set_visible(False)
            ax22.axis("off")

        label_pref = df.head(max_plot)[label]

    if out is not None:
        g.savefig(out)

    return fig
