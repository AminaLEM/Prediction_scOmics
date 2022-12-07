import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
from matplotlib.backends.backend_pdf import PdfPages
import sklearn.metrics as skm
import matplotlib.pyplot as plt
from collections import defaultdict

with open(snakemake.input[0], "rb") as b:
    test = pickle.load(b)

model_j_f = snakemake.input[1]
model_j_e = snakemake.input[2]
model_j_c = snakemake.input[3]
dim_exp = pd.read_csv(snakemake.input[5], index_col=0).shape[1] - 2
out1 = snakemake.output[0]
out2 = snakemake.output[1]
out3 = snakemake.output[2]
out_fig = snakemake.output[3]
out_j_txt = snakemake.output[4]
out_e_txt = snakemake.output[5]
out_c_txt = snakemake.output[6]


def confidence_interval(values):
    mean = np.mean(values)
    alpha = 0.95
    p = ((1.0 - alpha) / 2.0) * 100
    bottom = max(0.0, np.percentile(values, p))
    p = (alpha + ((1.0 - alpha) / 2.0)) * 100
    top = min(1.0, np.percentile(values, p))
    return mean, bottom, top


def eval_box(metrics, tit):
    f = plt.figure()
    val_vec = list()
    bt_vec = list()
    tp_vec = list()
    x_vec = list()
    for d in metrics:
        val, bt, tp = confidence_interval(metrics[d])
        val_vec.append(val)
        bt_vec.append(val - bt)
        tp_vec.append(tp - val)
        x_vec.append(d)

    with plt.rc_context({"figure.figsize": (4, 3), "figure.dpi": 300, "font.size": 16}):
        #         plt.errorbar(x_vec, val_vec, yerr=(bt_vec, tp_vec), linestyle="None",  fmt="ob",  capsize=3,  ecolor="k")
        plt.errorbar(
            x_vec,
            val_vec,
            yerr=(bt_vec, tp_vec),
            fmt="o",
            capsize=5,
            ecolor="k",
            lw=3,
            ls=":",
            color="blue",
        )
        plt.title(tit)
        plt.ylim(0.4, 1.05)
        plt.plot([], c="k", label="CI 95%")
        plt.plot([], c="blue", label="mean")
        plt.legend(loc="lower right")
    #         plt.savefig(tit, format='pdf', dpi=360)
    #         plt.show()
    return f


def compute_metrics(y):
    metrics = defaultdict(list)
    i = 0
    for l in y.columns:
        Ytest = test[i][:, -1]
        metrics["auc"].append(skm.roc_auc_score(Ytest, y[l]))
        metrics["acc"].append(skm.accuracy_score(Ytest, y[l] >= 0.5))
        metrics["f1"].append(skm.f1_score(Ytest, y[l] >= 0.5))
        metrics["rec"].append(skm.recall_score(Ytest, y[l] >= 0.5))
        metrics["prc"].append(skm.precision_score(Ytest, y[l] >= 0.5))
        metrics["auprc"].append(skm.average_precision_score(Ytest, y[l]))
        i = i + 1
    return metrics


if __name__ == "__main__":

    with open(model_j_f, "rb") as b:
        model_j = pickle.load(b)
    with open(model_j_e, "rb") as b:
        model_e = pickle.load(b)
    with open(model_j_c, "rb") as b:
        model_c = pickle.load(b)

    y_score1 = pd.DataFrame([])
    y_score2 = pd.DataFrame([])
    y_score3 = pd.DataFrame([])

    i = 0
    dim_cells = test[i].shape[1] - 4 - dim_exp
    for model in model_j:
        x_cell = test[i][:, dim_exp + 3 : (dim_exp + dim_cells + 3)]
        x_exp = test[i][:, 3 : dim_exp + 3]
        gender = test[i][:, :3]
        y_score1["sampling" + str(i)] = model.predict([x_exp, x_cell, gender]).flatten()
        i = i + 1
    i = 0
    for model in model_e:
        x_exp = test[i][:, : dim_exp + 3]
        y_score2["sampling" + str(i)] = model.predict(x_exp).flatten()
        i = i + 1
    i = 0
    for model in model_c:
        x_cell = np.concatenate(
            (test[i][:, :3], test[i][:, dim_exp + 3 : (dim_exp + 3 + dim_cells)]),
            axis=1,
        )
        y_score3["sampling" + str(i)] = model.predict(x_cell).flatten()
        i = i + 1

    y_score1.to_csv(out1)
    y_score2.to_csv(out2)
    y_score3.to_csv(out3)
    res_j = compute_metrics(y_score1)
    res_e = compute_metrics(y_score2)
    res_c = compute_metrics(y_score3)

    # compute mean and CI
    all_met = pd.DataFrame([])
    for d in res_j:
        val, bt, tp = confidence_interval(res_j[d])
        all_met[d] = [val, bt, tp]

    all_met.index = ["mean", "lower CI", "upper CI"]
    all_met.transpose().to_csv(out_j_txt)
    for d in res_e:
        val, bt, tp = confidence_interval(res_e[d])
        all_met[d] = [val, bt, tp]
    all_met.index = ["mean", "lower CI", "upper CI"]
    all_met.transpose().to_csv(out_e_txt)

    for d in res_c:
        val, bt, tp = confidence_interval(res_c[d])
        all_met[d] = [val, bt, tp]

    all_met.index = ["mean", "lower CI", "upper CI"]
    all_met.transpose().to_csv(out_c_txt)

    # plot figures
    fig1 = eval_box(
        {"CC": res_c["auc"], "GE": res_e["auc"], "CC&GE": res_j["auc"]}, "AUC"
    )
    fig2 = eval_box(
        {"CC": res_c["auprc"], "GE": res_e["auprc"], "CC&GE": res_j["auprc"]}, "AUPRC"
    )
    fig3 = eval_box(
        {"CC": res_c["acc"], "GE": res_e["acc"], "CC&GE": res_j["acc"]}, "ACCURACY"
    )
    fig4 = eval_box(
        {"CC": res_c["prc"], "GE": res_e["prc"], "CC&GE": res_j["prc"]}, "PRECISION"
    )
    fig5 = eval_box(
        {"CC": res_c["rec"], "GE": res_e["rec"], "CC&GE": res_j["rec"]}, "RECALL"
    )
    fig6 = eval_box(
        {"CC": res_c["f1"], "GE": res_e["f1"], "CC&GE": res_j["f1"]}, "F1-SCORE"
    )
    fig7 = eval_box(
        {"Precision": res_j["prc"], "Recall": res_j["rec"], "F1-score": res_j["f1"]},
        "Sens&Spec_CC&GE",
    )

    pp = PdfPages(out_fig)
    pp.savefig(fig1, bbox_inches="tight")
    pp.savefig(fig2, bbox_inches="tight")
    pp.savefig(fig3, bbox_inches="tight")
    pp.savefig(fig4, bbox_inches="tight")
    pp.savefig(fig5, bbox_inches="tight")
    pp.savefig(fig6, bbox_inches="tight")
    pp.savefig(fig7, bbox_inches="tight")

    pp.close()
