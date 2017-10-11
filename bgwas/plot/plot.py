"""
plotting function

:author: Benjamin Schubert
"""
import plotly
import seaborn as sns

def interactive_scatter(df, x_name, y_name, ann_name, x_title="", y_title=""):
    """
    generates an interactive scatter plot using plotly in ipython offline mode

    :param df: pandas dataframe which contains all relevant information
    :param x_name: the x column name
    :param y_name: the y column name
    :param ann_name: the annotation column name
    :param x_title: x axis title
    :param y_title: y axis title
    :return: plotly image
    """

    fig = {
        'data': [
            {
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