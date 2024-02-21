import pandas as pd
from optparse import OptionParser

####input paramater####
parser = OptionParser()
parser.add_option('-o', '--dir', help='work directory', dest='workdir', default='.')
(options, args) = parser.parse_args()

df = pd.read_csv(f"{options.workdir}/human_TE_cellline_all.csv", index_col=0)
transposed = df.T
transposed.index = transposed.index.str.replace("\.(.*)", "", regex=True)
transposed.to_csv(f"{options.workdir}/human_TE_cellline_all_T.csv")
