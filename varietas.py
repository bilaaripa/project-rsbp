import pandas as pd

class varietas:
    def __init__(self, sequence_path):
        self.data_frame = pd.read_csv(sequence_path)
        self.varietieslist = self.data_frame['N'].tolist()
        self.data_frame.drop('N', axis=1, inplace=True)

    def preprocessing_result(self):
        return self.data_frame