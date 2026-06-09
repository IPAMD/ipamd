"""plot the data"""
from functools import singledispatch
from matplotlib import pyplot as plt
import numpy as np

from ipamd.public.models.data import PointSet, Vector, Matrix, Ratio, Distribution, Scalar, String
from ipamd.public.utils.output import warning
from ipamd.public.utils.plugin_manager_v1 import PluginBase


@singledispatch
def plot(data, **kwargs):
    """
    plot the data
    
    :param data: data to plot
    :param kwargs: keyword arguments
    :return: None
    """
    raise NotImplementedError()

@plot.register
def _(data: PointSet, **kwargs):
    x_data = [point[0] for point in data.data]
    y_data = [point[1] for point in data.data]
    plt.scatter(x_data, y_data)
    plt.xlabel(data.meta['x_label'])
    plt.ylabel(data.meta['y_label'])

@plot.register
def _(data: Vector, **kwargs):
    size_x = len(data.data)

    appearance = kwargs.get('appearance', 'line')
    if appearance == 'bar':
        plt.bar(data.meta['x_axis'], data.data)
    elif appearance == 'heatmap':
        plt.imshow([data.data], aspect='auto')
        plt.yticks([])
        plt.colorbar()
    elif appearance == 'discrete_heatmap':
        unique_categories = np.unique(data.data)
        n_categories = len(unique_categories)

        category_to_num = {cat: i for i, cat in enumerate(unique_categories)}
        numeric_data = np.vectorize(category_to_num.get)(data.data)
        discrete_cmap = plt.cm.get_cmap('viridis', n_categories)
        im = plt.imshow(
            [numeric_data],
            aspect='auto',
            cmap=discrete_cmap,
            vmin=-0.5,
            vmax=n_categories - 0.5
        )
        cbar = plt.colorbar(im, ticks=range(n_categories))
        cbar.set_ticklabels(unique_categories)  # 将数值标签映射回原始类别文本
        cbar.set_label(data.meta.get('z_label', 'Categories'))
        plt.yticks([])
    else:
        plt.plot(data.meta['x_axis'], data.data)
    is_series = kwargs.get('is_series', False)
    if is_series:
        label_interval_x = max(1, size_x // 10)
        plt.xticks(
            ticks=list(range(0, size_x, label_interval_x)),
            labels=data.meta['x_axis'][::label_interval_x]
        )

    plt.xlabel(data.meta['x_label'])
    plt.ylabel(data.meta['y_label'])

@plot.register
def _(data: Ratio, **kwargs):
    plt.pie(
        data.data,
        labels=[
            label
            if value / sum(data.data) > 0.02
            else f'{label} - {value / sum(data.data) * 100:.1f}%'
            for value, label
            in zip(data.data, data.meta['labels'])
        ],
        autopct=lambda pct: f'{pct:.1f}%' if pct > 2 else '',
        startangle=90,
        shadow=True,
        pctdistance=0.75,
        explode=[0.1 if i < max(data.data) else 0 for i in data.data]
    )

@plot.register
def _(data: Distribution, **kwargs):
    plt.hist(data.data, data.meta['bins'], edgecolor='black')

@plot.register
def _(data: Matrix, **kwargs):
    plt.imshow(
        data.data,
        vmin=kwargs.get('vmin', None),
        vmax=kwargs.get('vmax', None),
        origin='lower',
        aspect='auto'
    )
    plt.colorbar()
    size_x = data.data.shape[1]
    size_y = data.data.shape[0]
    label_interval_x = max(1, size_x // 10)
    label_interval_y = max(1, size_y // 10)
    plt.xlabel(data.meta['x_label'])
    plt.xticks(
        ticks=list(range(0, size_x, label_interval_x)),
        labels=data.meta['x_axis'][::label_interval_x]
    )
    plt.ylabel(data.meta['y_label'])
    plt.yticks(
        ticks=list(range(0, size_y, label_interval_y)),
        labels=data.meta['y_axis'][::label_interval_y]
    )

def func(data, style=None, save_figure=False, **kwargs):
    """plugin main function
    :param data: data to plot
    :param style: style of the plot
    :param save_figure: whether to save the figure
    :param kwargs: keyword arguments
    :return: None
    """
    if isinstance(data, Scalar) or isinstance(data, String):
        warning("data cannot be plotted.")
        PluginBase.call('print', data)
        return

    with plt.rc_context(style):

        plot(data, **kwargs)
        plt.title(data.meta['title'])
        plt.tight_layout()

        if save_figure:
            plt.savefig(f"{data.meta['title']}.png")
        plt.show()
