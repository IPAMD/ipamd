import numpy as np

class AnalysisResult:
    def __init__(self, title, data):
        self.meta = {
            "title": title,
        }
        self.data = np.array(data)

    def update(self, new):
        self.data = np.array(new)

class Scalar(AnalysisResult):
    def __init__(self, title, data, unit=""):
        super().__init__(title, data)
        self.meta.update({
            "unit": unit,
        })

class Vector(AnalysisResult):
    def __init__(self, title, data, x_label="Keys", y_label='Values', x_axis=None):
        super().__init__(title, data)
        self.meta.update({
            "x_label": x_label,
            "y_label": y_label,
            "x_axis": x_axis if x_axis is not None else list(map(str, list(range(len(self.data)))))
        })

class PointSet(AnalysisResult):
    def __init__(self, title, data, x_label="x", y_label="y"):
        super().__init__(title, data)

        self.meta.update({
            "x_label": x_label,
            "y_label": y_label
        })

class Distribution(AnalysisResult):
    def __init__(self, title, data, bins=10, label="Frequency"):
        super().__init__(title, data)

        self.meta.update({
            "bins": bins,
            "label": label
        })

class Matrix(AnalysisResult):
    def __init__(self, title, data, x_label="Columns", y_label="Rows"):
        super().__init__(title, data)

        self.meta.update({
            "x_label": x_label,
            "y_label": y_label,
            "x_axis": list(map(str, list(range(self.data.shape[1])))),
            "y_axis": list(map(str, list(range(self.data.shape[0])))),
        })

class Ratio(AnalysisResult):
    def __init__(self, title, data, labels):
        super().__init__(title, data)

        self.meta.update({
            "labels": labels
        })

class String(AnalysisResult):
    """
    String data type
    """
    def __init__(self, title, data):
        self.meta = {
            "title": title
        }
        self.data = data
