import numpy as np

class Matfile:

    def __init__(self, path):
        self.rows = -1
        self.columns = -1

        header_lines = []
        data_lines = []
        with open(path, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    header_lines.append(line)
                else:
                    data_lines.append(line)
        self.read_header(header_lines)
        self.validate_header()
        self.read_data(data_lines)

    def read_header(self, lines):
        for line in lines:
            (key, value) = [tok.strip() for tok in line.strip("# ").split(":")]
            if key.lower() == "name":
                self.name = value
            elif key.lower() == "rows":
                self.rows = int(value)
            elif key.lower() == "columns":
                self.columns = int(value)
            elif key.lower() == "type":
                if value.lower() != "matrix":
                    raise Exception("Only matrix type is supported")

    def validate_header(self):
        if self.rows < 1:
            raise Exception("Header has missing or invalid 'rows' field")
        elif self.columns < 1:
            raise Exception("Header has missing or invalid 'columns' field")

    def read_data(self, lines):
        self.array = np.empty((self.rows, self.columns), np.float32)
        for i, line in enumerate(lines):
            ents = [float(tok) for tok in line.split()]
            if len(ents) != self.columns:
                raise Exception("Bad line length")
            for j in range(len(ents)):
                self.array[i][j] = ents[j]

    def get_ndarray(self):
        return self.array
