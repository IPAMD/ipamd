from functools import singledispatch
import numpy as np
from ipamd.public.models.data import *
from ipamd.public.utils.output import output, tabulate


configure = {
    "alias": 'print',
}

def func(data, precision=3):
    output(f'{data.meta["title"]}:')
    print_data(data, precision)

@singledispatch
def print_data(data, precision):
    raise NotImplementedError(f'Output for data type {type(data)} is not implemented.')

@print_data.register
def _(data: Scalar, precision):
    output(f'{data.data:.{precision}f} {data.meta['unit']}')

@print_data.register
def _(data: Vector, precision):
    points = []
    for i in range(len(data.data)):
        x = data.meta['axis'][i]
        y = data.data[i]
        points.append(f'({x:.{precision}f}, {y:.{precision}f})')
    output(' '.join(points))

@print_data.register
def _(data: Vector, precision):
    tabulate('', data.meta['x_axis'], [[f'{y:.{precision}f}' for y in data.data]])

@print_data.register
def _(data: PointSet, precision):
    points = []
    for point in data.data:
        points.append([f'{point[0]:.{precision}f}', f'{point[1]:.{precision}f}'])
    tabulate('', [data.meta['x_label'], data.meta['y_label']], points)

@print_data.register
def _(data: Matrix, precision):
    rows = []
    for row, label in zip(data.data, data.meta['y_axis']):
        rows.append([label, *[f'{x:.{precision}f}' for x in row]])
    tabulate('', ['*', *data.meta['x_axis']], rows)

@print_data.register
def _(data: Ratio, precision):
    tabulate('', data.meta['labels'], [[f'{x:.{precision}f}' for x in data.data]])

@print_data.register
def _(data: Distribution, precision):
    _, edges = np.histogram(data.data, bins=data.meta['bins'])
    for i in range(data.meta['bins']):
        if i == data.meta['bins'] - 1:
            bin_range = f'[{edges[i]:.{precision}f}, {edges[i+1]:.{precision}f}]'
            group_i = data.data[(data.data >= edges[i]) & (data.data <= edges[i + 1])]
        else:
            bin_range = f'[{edges[i]:.{precision}f}, {edges[i+1]:.{precision}f})'
            group_i = data.data[(data.data >= edges[i]) & (data.data < edges[i + 1])]

        formatted_row = ', '.join([f'{x:.{precision}f}' for x in group_i])
        output(f'{bin_range}: [{formatted_row}]')

@print_data.register
def _(data: String, precision):
    output(f'{data.meta["title"]}: "{data.data}"')
