
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as ticker


def get_plot(pvl, cdf, k , type):
    plt.plot(pvl, cdf, marker='.', linestyle='none')
    plt.xlim(0, 0.01)
    plt.ylim(0, 0.2)
    plt.xlabel('P-values')
    # plt.set_xlim(x)
    plt.ylabel('Cumulative Probability')
    plt.title('Cumulative Distribution Function (CDF) of P-values')
    # plt.grid(True)
    plt.savefig(f'{k}_{type}_plot.svg',format='svg')

  
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd



# def get_plot_multiple(data_sets, labels, k, plot_type):
#     plt.figure(figsize=(8, 6))  
#     for i, data in enumerate(data_sets):
#         plt.plot(data['p-Value'], data['CDF'], marker='.', linestyle='none', label=labels[i])
    
#     plt.xlim(0, 0.01)
#     plt.ylim(0, 0.2)
#     plt.xlabel('P-values')
#     plt.ylabel('Cumulative Probability')
#     plt.title('Cumulative Distribution Function (CDF) of P-values')
#     plt.legend()
#     plt.savefig(f'{k}_{plot_type}_plot.svg', format='svg')

# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd

# def get_plot_multiple(data_sets, labels, k, plot_type):
#     plt.figure(figsize=(8, 6))  
#     for i, data in enumerate(data_sets):
#         plt.plot(data['p-Value'], data['CDF'], label=labels[i], linewidth=2.5)
    
#     plt.xlim(0, 0.01)
#     plt.ylim(0, 0.2)
#     plt.xlabel('P-values')
#     plt.ylabel('Cumulative Probability')
#     plt.title('Cumulative Distribution Function (CDF) of P-values')
#     plt.legend()
#     plt.savefig(f'{k}_{plot_type}_plot.svg', format='svg')


# get_plot_multiple([data1, data2], ['Label1', 'Label2'], 'K', 'Type')
