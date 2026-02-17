import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd



def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

# 读取数据

s_file_path = r"D:\data\lstSimulate\Result\\"

methods = ["RF","LSTM"]   #,"LSTM"
indicators = ["RMSE","r2","bias"]  #,"MAE","r2","bias"
stations = ["AR","DM","SDQ"]

for method in methods:
    for indicator in indicators:
        file_path = r"D:\data\lstSimulate\Result/" + method + r"\indicator/"
        save_file_path = s_file_path + method + r"\graph\\" + indicator + r"\\"

        # # 1：十个普通站
        # data_path= file_path + "1/" + indicator + ".xlsx"
        # data = pd.read_excel(data_path,sheet_name='1',header=0,index_col=0)
        # row_label = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'ALL']
        # col_label = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'ALL']
        # fig, ax = plt.subplots(figsize=(12, 9))
        # ax.set_title(method + " regression "+ indicator)
        # plt.xlabel("Train Station")
        # plt.ylabel("Test Station")
        # threshold = 0
        # if indicator == "RMSE":
        #     vmin,vmax = 0,6
        # elif indicator == "MAE":
        #     vmin, vmax = -5,5
        # elif indicator == "r2":
        #     vmin, vmax = 0.80,1
        # elif indicator == 'bias':
        #     vmin , vmax , threshold = -5,5,-10
        # im, cbar = heatmap(data, row_label, col_label, ax=ax,
        #                    cmap='rainbow', cbarlabel=indicator, vmin=vmin, vmax=vmax)
        # texts = annotate_heatmap(im, valfmt="{x:.2f} ", textcolors=("white", "black"), threshold=threshold)
        # fig.tight_layout()
        # path = save_file_path +"1_" +  method + "_ "+ indicator + ".png"
        # plt.savefig(path, dpi=300)
        # # plt.show()
        #
        # # 2：不同植被类型
        # data_path= file_path + "2/" + indicator + ".xlsx"
        # data = pd.read_excel(data_path,sheet_name='1',header=0,index_col=0)
        # col_label = ['M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'ALL']
        # row_label = ['M11', 'M12', 'M13','ALL']
        # fig, ax = plt.subplots(figsize=(12,4))
        # ax.set_title(method + " regression "+ indicator)
        # plt.xlabel("Train Station")
        # plt.ylabel("Test Station")
        # threshold = 0
        # if indicator == "RMSE":
        #     vmin,vmax = 0,6
        # elif indicator == "MAE":
        #     vmin, vmax = -5,5
        # elif indicator == "r2":
        #     vmin, vmax = 0.80,1
        # elif indicator == 'bias':
        #     vmin , vmax , threshold = -5,5,-10
        # im, cbar = heatmap(data, row_label, col_label, ax=ax,
        #                    cmap='rainbow', cbarlabel=indicator, vmin=vmin, vmax=vmax)
        # texts = annotate_heatmap(im, valfmt="{x:.2f} ", textcolors=("white", "black"), threshold=threshold)
        # fig.tight_layout()
        # path = save_file_path +"2_" +  method + "_" + indicator + ".png"
        # plt.savefig(path, dpi=300)
        # # plt.show()
        # #
        # # 3：同站点，不同年
        # for station in stations:
        #     data_path = file_path + "3/3_"+station+"/" + indicator + ".xlsx"
        #     data = pd.read_excel(data_path, sheet_name='1', header=0, index_col=0)
        #     col_label = [station + '\n2013', station +'\n2014', station +'\n2015', station +'\n2016', station +'\n2017', 'ALL']
        #     row_label = [station + '\n2013', station +'\n2014', station +'\n2015', station +'\n2016', station +'\n2017', 'ALL']
        #     fig, ax = plt.subplots(figsize=(6,5))
        #     ax.set_title(method + " regression " + indicator)
        #     plt.xlabel("Train Station")
        #     plt.ylabel("Test Station")
        #     threshold = 0
        #     if indicator == "RMSE":
        #         vmin, vmax = 0, 6
        #     elif indicator == "MAE":
        #         vmin, vmax = -5, 5
        #     elif indicator == "r2":
        #         vmin, vmax = 0.80, 1
        #     elif indicator == 'bias':
        #         vmin, vmax, threshold = -5, 5, -10
        #     im, cbar = heatmap(data, row_label, col_label, ax=ax,
        #                        cmap='rainbow', cbarlabel=indicator, vmin=vmin, vmax=vmax)
        #     texts = annotate_heatmap(im, valfmt="{x:.2f} ", textcolors=("white", "black"), threshold=threshold)
        #     fig.tight_layout()
        #     path = save_file_path + "3_" + method + "_" +station + "_"+ indicator + ".png"
        #     plt.savefig(path, dpi=300)
        #     # plt.show()


        # # 4：不同站点，不同年
        # trainstation = ["DM"]
        # teststation = ["SDQ","AR"]
        # for tstation in teststation:
        #     data_path = file_path + "4/4_" + trainstation[0] + "_"+tstation + "/" + indicator + ".xlsx"
        #     data = pd.read_excel(data_path, sheet_name='1', header=0, index_col=0)
        #     row_label = [tstation + '\n2013', tstation + '\n2014', tstation + '\n2015', tstation + '\n2016',
        #                  tstation + '\n2017', 'ALL']
        #     col_label = [trainstation[0] + '\n2013', trainstation[0] + '\n2014', trainstation[0] + '\n2015', trainstation[0] + '\n2016',
        #                  trainstation[0] + '\n2017', 'ALL']
        #     fig, ax = plt.subplots(figsize=(6, 5))
        #     ax.set_title(method + " regression " + indicator)
        #     plt.xlabel("Train Station")
        #     plt.ylabel("Test Station")
        #     threshold = 0
        #     if indicator == "RMSE":
        #         vmin, vmax = 0, 6
        #     elif indicator == "MAE":
        #         vmin, vmax = -5, 5
        #     elif indicator == "r2":
        #         vmin, vmax = 0.80, 1
        #     elif indicator == 'bias':
        #         vmin, vmax, threshold = -5, 5, -10
        #     im, cbar = heatmap(data, row_label, col_label, ax=ax,
        #                        cmap='rainbow', cbarlabel=indicator, vmin=vmin, vmax=vmax)
        #     texts = annotate_heatmap(im, valfmt="{x:.2f} ", textcolors=("white", "black"), threshold=threshold)
        #     fig.tight_layout()
        #     path = save_file_path + "4_" + method + "_" + trainstation[0] + "_" +tstation+"_"+ indicator + ".png"
        #     plt.savefig(path, dpi=300)
        #
        #     # plt.show()

        # 5.4的总结图
        data_path = "D:\data\lstSimulate\Result\\" + method +r"\indicator\4\\" + indicator + ".xlsx"
        for i in [1,2,3]:
            if i == 1:
                sname,tname1, tname2 = "DM","AR","SDQ"
            elif i ==2:
                sname, tname1, tname2 = "AR","DM" , "SDQ"
            else:
                sname, tname1, tname2 = "SDQ", "DM", "AR"
            data = pd.read_excel(data_path, sheet_name= str(i), header=0, index_col=0)
            col_label = [sname + '\n2013', sname + '\n2014', sname + '\n2015', sname + '\n2016', sname + '\n2017', sname + '\n2013-17']
            row_label = [tname1,tname2]
            fig, ax = plt.subplots(figsize=(6, 2))
            ax.set_title(method + " regression " + indicator)
            plt.xlabel("Test Station")
            plt.ylabel("Train Station")
            threshold = 0
            if indicator == "RMSE":
                vmin, vmax = 0, 6
            elif indicator == "MAE":
                vmin, vmax = -5, 5
            elif indicator == "r2":
                vmin, vmax = 0.80, 1
            elif indicator == 'bias':
                vmin, vmax, threshold = -5, 5, -10
            im, cbar = heatmap(data, row_label, col_label, ax=ax,
                               cmap='rainbow', cbarlabel=indicator, vmin=vmin, vmax=vmax)
            texts = annotate_heatmap(im, valfmt="{x:.2f} ", textcolors=("white", "black"), threshold=threshold)
            fig.tight_layout()
            path = save_file_path + "5_" + method + "_" + indicator +"_" + sname+ ".png"
            plt.savefig(path, dpi=300)
            # plt.show()

        #6
        data_path = "D:\data\lstSimulate\Result\\" + method + r"\indicator\4\\"+indicator+".xlsx"
        data = pd.read_excel(data_path,sheet_name='4',header=0,index_col=0)
        col_label = ['DM', 'AR', 'SDQ']
        row_label = ['DM', 'AR', 'SDQ']
        fig, ax = plt.subplots(figsize=(4,4))
        ax.set_title(method + " regression "+ indicator)
        plt.xlabel("Train Station")
        plt.ylabel("Test Station")
        threshold = 0
        if indicator == "RMSE":
            vmin,vmax = 0,6
        elif indicator == "MAE":
            vmin, vmax = -5,5
        elif indicator == "r2":
            vmin, vmax = 0.80,1
        elif indicator == 'bias':
            vmin , vmax , threshold = -5,5,-10
        im, cbar = heatmap(data, row_label, col_label, ax=ax,
                           cmap='rainbow', cbarlabel=indicator, vmin=vmin, vmax=vmax)
        texts = annotate_heatmap(im, valfmt="{x:.2f} ", textcolors=("white", "black"), threshold=threshold)
        fig.tight_layout()
        path = save_file_path +"6_" +  method + "_" + indicator + ".png"
        plt.savefig(path, dpi=300)
        # plt.show()



