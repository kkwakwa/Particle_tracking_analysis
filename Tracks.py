import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class Tracks:

    def __init__(self, filepath):
        self.trackdata = pd.read_csv(filepath)