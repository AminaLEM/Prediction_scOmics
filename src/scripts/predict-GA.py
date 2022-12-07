from sklearn.preprocessing import LabelEncoder, minmax_scale
import pickle
import numpy as np
import pandas as pd
import tensorflow as tf
from matplotlib.backends.backend_pdf import PdfPages
import sklearn.metrics as skm
import matplotlib.pyplot as plt
from collections import defaultdict
import shap
from numba import njit, prange
import seaborn as sns
from statannot import add_stat_annotation

# For SHAP package should disable v2 tensorflow
tf.compat.v1.disable_v2_behavior()
x_cell = pd.DataFrame([])
for elem in snakemake.input["CC"]:
    x_ = pd.read_csv(elem, index_col=0)
    x_cell = pd.concat([x_cell, x_], axis=0)

x_exp = pd.DataFrame([])
for elem in snakemake.input["GE"]:
    x_ = pd.read_csv(elem, index_col=0)
    x_exp = pd.concat([x_exp, x_], axis=0)

# Load data and models
model_j_f = snakemake.input["model_j"]
model_j_e = snakemake.input["model_e"]
model_j_c = snakemake.input["model_c"]
svm_e_f = snakemake.input["svm_e"]
svm_c_f = snakemake.input["svm_c"]
LogReg_e_f = snakemake.input["LogReg_e"]
LogReg_c_f = snakemake.input["LogReg_c"]
RF_e_f = snakemake.input["RF_e"]
RF_c_f = snakemake.input["RF_c"]
svm_j_f = snakemake.input["svm_j"]
LogReg_j_f = snakemake.input["LogReg_j"]
RF_j_f = snakemake.input["RF_j"]

out_fig = snakemake.output["fig"]
svm_fig = snakemake.output["svm_fig"]
LogReg_fig = snakemake.output["LogReg_fig"]
RF_fig = snakemake.output["RF_fig"]
out_shap = snakemake.output["out_shap"]


gender = pd.read_csv(snakemake.input["gender"], index_col=0)
gender_dummies = pd.get_dummies(gender.gender, prefix="gender_")
gender_dummies.index = gender.index
age_ = pd.read_csv(snakemake.input["age"], index_col=0)
age_.age = age_.age.astype("category").cat.codes / len(age_.age.unique())


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

    return f


def eval_box_sbn(metrics, tit):
    f = plt.figure()
    data = pd.DataFrame([], columns=["modality", "value"])
    for d in metrics:
        data_ = pd.DataFrame(
            [np.full((len(metrics[d])), d), metrics[d]], index=["modality", "value"]
        )
        data = pd.concat([data, data_.transpose()], ignore_index=True)
    data.value = data.value.astype("float")
    mod = data.modality.unique()
    if len(mod) == 3:
        pairs = [("CC", "GE"), ("CC", "CC&GE"), ("GE", "CC&GE")]
    else:
        pairs = [
            ("CC&GE", "SVM"),
            ("CC&GE", "RF"),
            ("CC&GE", "LogisticR"),
        ]
    with plt.rc_context({"figure.figsize": (10, 8), "figure.dpi": 96, "font.size": 16}):

        g = sns.boxplot(x="modality", y="value", data=data)

        add_stat_annotation(
            g,
            x="modality",
            y="value",
            data=data,
            box_pairs=pairs,
            test="Mann-Whitney",
            text_format="star",
            loc="inside",
            verbose=2,
        )
        plt.ylabel(tit)
        if len(mod) == 3:
            plt.xlabel("modality")
        else:
            plt.xlabel("model")

    return f


def compute_metrics(y):
    metrics = defaultdict(list)
    col_len = y.shape[1]
    for i in range(col_len):
        l = y.columns[i]
        metrics["auc"].append(skm.roc_auc_score(Ytest, y[l]))
        metrics["acc"].append(skm.accuracy_score(Ytest, y[l] >= 0.5))
        metrics["f1"].append(skm.f1_score(Ytest, y[l] >= 0.5))
        metrics["rec"].append(skm.recall_score(Ytest, y[l] >= 0.5))
        metrics["prc"].append(skm.precision_score(Ytest, y[l] >= 0.5))
        metrics["auprc"].append(skm.average_precision_score(Ytest, y[l]))
    return metrics


def predict_loop(model_, data):
    y_score1 = pd.DataFrame([])
    model_len = len(model_)
    for i in range(model_len):
        model = model_[i]
        y_score1["sampling" + str(i)] = model.predict(data).flatten()
    return y_score1


def shap_loop(model_j, training_set, dim_exp, dim_cells, x_exp, x_cell):
    model_len = len(model_j)
    for i in range(model_len):
        model_joint = model_j[i]
        x_ref_exp = training_set[i][:, 3 : dim_exp + 3]
        x_ref_cell = training_set[i][:, dim_exp + 3 : (dim_exp + dim_cells + 3)]
        x_ref_gender = training_set[i][:, :3]
        explainer = shap.DeepExplainer(
            model_joint, [x_ref_exp, x_ref_cell, x_ref_gender]
        )
        shap_values = explainer.shap_values(
            [
                np.array(x_exp),
                np.array(x_cell),
                np.concatenate((gender_dummies, age_), axis=1),
            ]
        )
        if model_joint == model_j[0]:
            shap_values_all_exp = shap_values[0][0]
            shap_values_all_cell = shap_values[0][1]
            shap_values_all_gender = shap_values[0][2]
        else:
            shap_values_all_exp = shap_values_all_exp + shap_values[0][0]
            shap_values_all_cell = shap_values_all_cell + shap_values[0][1]
            shap_values_all_gender = shap_values_all_gender + shap_values[0][2]
    return shap_values_all_exp, shap_values_all_cell, shap_values_all_gender


if __name__ == "__main__":

    with open(model_j_f, "rb") as b:
        model_j = pickle.load(b)
    with open(model_j_e, "rb") as b:
        model_e = pickle.load(b)
    with open(model_j_c, "rb") as b:
        model_c = pickle.load(b)
    with open(svm_j_f, "rb") as b:
        svm_j = pickle.load(b)
    with open(svm_e_f, "rb") as b:
        svm_e = pickle.load(b)
    with open(svm_c_f, "rb") as b:
        svm_c = pickle.load(b)

    with open(LogReg_j_f, "rb") as b:
        LogReg_j = pickle.load(b)
    with open(LogReg_e_f, "rb") as b:
        LogReg_e = pickle.load(b)
    with open(LogReg_c_f, "rb") as b:
        LogReg_c = pickle.load(b)
    with open(RF_j_f, "rb") as b:
        RF_j = pickle.load(b)
        
    with open(RF_e_f, "rb") as b:
        RF_e = pickle.load(b)
    with open(RF_c_f, "rb") as b:
        RF_c = pickle.load(b)

    # prepare data
    x_exp = x_exp.loc[x_exp["condition"].isin(["Mild", "Severe"]), :]
    x_cell = x_cell.loc[x_cell["condition"].isin(["Mild", "Severe"]), :]
    x_cell = x_cell.loc[x_exp.index, :]
    gender_dummies = gender_dummies.loc[x_exp.index.astype(str), :]
    age_ = age_.loc[x_exp.index.astype(str), :]
    label = x_cell.iloc[:, -1].values
    x_cell = x_cell.drop("condition", axis=1)
    x_exp = x_exp.drop("condition", axis=1)
    x_exp = x_exp.drop("who_score", axis=1)
    genes = x_exp.columns
    cells = x_cell.columns
    selected_cols = cells

    le = LabelEncoder()
    Ytest = le.fit_transform(label)
    x_exp = x_exp / 15
    x_cell = x_cell.div(x_cell.sum(axis=1), axis=0)
    x_cell = x_cell.loc[:, selected_cols]

    # run prediction
    y_score1_RF = predict_loop(
        RF_j, np.concatenate((gender_dummies, age_, x_exp, x_cell), axis=1)
    )
    y_score2_RF = predict_loop(
        RF_e, np.concatenate((gender_dummies, age_, x_exp), axis=1)
    )
    y_score3_RF = predict_loop(
        RF_c, np.concatenate((gender_dummies, age_, x_cell), axis=1)
    )

    y_score1 = predict_loop(
        model_j, [x_exp, x_cell, np.concatenate((gender_dummies, age_), axis=1)]
    )
    y_score2 = predict_loop(
        model_e, np.concatenate((gender_dummies, age_, x_exp), axis=1)
    )
    y_score3 = predict_loop(
        model_c, np.concatenate((gender_dummies, age_, x_cell), axis=1)
    )

    y_score1_svm = predict_loop(
        svm_j, np.concatenate((gender_dummies, age_, x_exp, x_cell), axis=1)
    )
    y_score2_svm = predict_loop(
        svm_e, np.concatenate((gender_dummies, age_, x_exp), axis=1)
    )
    y_score3_svm = predict_loop(
        svm_c, np.concatenate((gender_dummies, age_, x_cell), axis=1)
    )

    y_score1_LogReg = predict_loop(
        LogReg_j, np.concatenate((gender_dummies, age_, x_exp, x_cell), axis=1)
    )
    y_score2_LogReg = predict_loop(
        LogReg_e, np.concatenate((gender_dummies, age_, x_exp), axis=1)
    )
    y_score3_LogReg = predict_loop(
        LogReg_c, np.concatenate((gender_dummies, age_, x_cell), axis=1)
    )

    # compute metrics and plot results
    all_yscores = {
        out_fig: [y_score1, y_score2, y_score3],
        svm_fig: [y_score1_svm, y_score2_svm, y_score3_svm],
        LogReg_fig: [y_score1_LogReg, y_score2_LogReg, y_score3_LogReg],
        RF_fig: [y_score1_RF, y_score2_RF, y_score3_RF],
    }
    comp = {}
    for key in all_yscores.keys():
        all_yscores[key][0].to_csv(key.replace(".pdf", "_CC_GE.csv"))
        all_yscores[key][1].to_csv(key.replace(".pdf", "_GE.csv"))
        all_yscores[key][2].to_csv(key.replace(".pdf", "_CC.csv"))
        res_j = compute_metrics(all_yscores[key][0])
        res_e = compute_metrics(all_yscores[key][1])
        res_c = compute_metrics(all_yscores[key][2])
        comp[key] = res_j
        # compute mean and CI
        all_met = pd.DataFrame([])
        for d in res_j:
            val, bt, tp = confidence_interval(res_j[d])
            all_met[d] = [val, bt, tp]

        all_met.index = ["mean", "lower CI", "upper CI"]
        all_met.transpose().to_csv(key.replace(".pdf", "_CC_GE.txt"))
        for d in res_e:
            val, bt, tp = confidence_interval(res_e[d])
            all_met[d] = [val, bt, tp]
        all_met.index = ["mean", "lower CI", "upper CI"]
        all_met.transpose().to_csv(key.replace(".pdf", "_GE.txt"))

        for d in res_c:
            val, bt, tp = confidence_interval(res_c[d])
            all_met[d] = [val, bt, tp]

        all_met.index = ["mean", "lower CI", "upper CI"]
        all_met.transpose().to_csv(key.replace(".pdf", "_CC.txt"))

        # plot figures
        fig1 = eval_box(
            {"CC": res_c["auc"], "GE": res_e["auc"], "CC&GE": res_j["auc"]}, "AUC"
        )
        fig2 = eval_box(
            {"CC": res_c["auprc"], "GE": res_e["auprc"], "CC&GE": res_j["auprc"]},
            "AUPRC",
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
            {
                "Precision": res_j["prc"],
                "Recall": res_j["rec"],
                "F1-score": res_j["f1"],
            },
            "Sens&Spec_CC&GE",
        )

        fig8 = eval_box_sbn(
            {"CC": res_c["auc"], "GE": res_e["auc"], "CC&GE": res_j["auc"]}, "AUC"
        )
        fig9 = eval_box_sbn(
            {"CC": res_c["auprc"], "GE": res_e["auprc"], "CC&GE": res_j["auprc"]},
            "Average Precision",
        )
        fig10 = eval_box_sbn(
            {"CC": res_c["acc"], "GE": res_e["acc"], "CC&GE": res_j["acc"]}, "Accuracy"
        )
        fig11 = eval_box_sbn(
            {"CC": res_c["prc"], "GE": res_e["prc"], "CC&GE": res_j["prc"]}, "Precision"
        )
        fig12 = eval_box_sbn(
            {"CC": res_c["rec"], "GE": res_e["rec"], "CC&GE": res_j["rec"]}, "Recall"
        )
        fig13 = eval_box_sbn(
            {"CC": res_c["f1"], "GE": res_e["f1"], "CC&GE": res_j["f1"]}, "F1-Score"
        )

        pp = PdfPages(key)
        pp.savefig(fig1, bbox_inches="tight")
        pp.savefig(fig2, bbox_inches="tight")
        pp.savefig(fig3, bbox_inches="tight")
        pp.savefig(fig4, bbox_inches="tight")
        pp.savefig(fig5, bbox_inches="tight")
        pp.savefig(fig6, bbox_inches="tight")
        pp.savefig(fig7, bbox_inches="tight")
        pp.savefig(fig8, bbox_inches="tight")
        pp.savefig(fig9, bbox_inches="tight")
        pp.savefig(fig10, bbox_inches="tight")
        pp.savefig(fig11, bbox_inches="tight")
        pp.savefig(fig12, bbox_inches="tight")
        pp.savefig(fig13, bbox_inches="tight")
        pp.close()

    fig14 = eval_box_sbn(
        {
            "LogisticR": comp[LogReg_fig]["auc"],
            "SVM": comp[svm_fig]["auc"],
            "RF": comp[RF_fig]["auc"],
            "CC&GE": comp[out_fig]["auc"],
        },
        "AUC",
    )
    fig15 = eval_box_sbn(
        {
            "LogisticR": comp[LogReg_fig]["auprc"],
            "SVM": comp[svm_fig]["auprc"],
            "RF": comp[RF_fig]["auprc"],
            "CC&GE": comp[out_fig]["auprc"],
        },
        "Average Precision",
    )
    fig16 = eval_box_sbn(
        {
            "LogisticR": comp[LogReg_fig]["acc"],
            "SVM": comp[svm_fig]["acc"],
            "RF": comp[RF_fig]["acc"],
            "CC&GE": comp[out_fig]["acc"],
        },
        "Accuracy",
    )
    fig17 = eval_box_sbn(
        {
            "LogisticR": comp[LogReg_fig]["prc"],
            "SVM": comp[svm_fig]["prc"],
            "RF": comp[RF_fig]["prc"],
            "CC&GE": comp[out_fig]["prc"],
        },
        "Precision",
    )
    fig18 = eval_box_sbn(
        {
            "LogisticR": comp[LogReg_fig]["rec"],
            "SVM": comp[svm_fig]["rec"],
            "RF": comp[RF_fig]["rec"],
            "CC&GE": comp[out_fig]["rec"],
        },
        "Recall",
    )
    fig19 = eval_box_sbn(
        {
            "LogisticR": comp[LogReg_fig]["f1"],
            "SVM": comp[svm_fig]["f1"],
            "RF": comp[RF_fig]["f1"],
            "CC&GE": comp[out_fig]["f1"],
        },
        "F1-Score",
    )

    all_yscores_ = defaultdict(list)
    for key in all_yscores.keys():
        all_yscores_[key.split("/")[-1] + "_GE_CC"] = all_yscores[key][0].mean(axis=1)
        if key == out_fig:
            all_yscores_[key.split("/")[-1] + "GE"] = all_yscores[key][1].mean(axis=1)
            all_yscores_[key.split("/")[-1] + "CC"] = all_yscores[key][2].mean(axis=1)

    fig20 = plt.figure()
    with plt.rc_context({"figure.figsize": (6, 5), "figure.dpi": 300, "font.size": 10}):
        for key in all_yscores_.keys():
            ns_auc = skm.roc_auc_score(Ytest, all_yscores_[key])
            # calculate roc curves
            fpr, tpr, _ = skm.roc_curve(Ytest, all_yscores_[key])
            # plot the roc curve for the model
            plt.plot(fpr, tpr, marker=".", label=key + ": " + str(round(ns_auc, 2)))
            # axis labels
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            # show the legend
            plt.legend()
            # show the plot
    #         plt.show()

    fig21 = plt.figure()
    with plt.rc_context({"figure.figsize": (6, 5), "figure.dpi": 300, "font.size": 10}):
        for key in all_yscores_.keys():

            precision, recall, _ = skm.precision_recall_curve(Ytest, all_yscores_[key])
            auc = skm.auc(recall, precision)
            avr_prec = skm.average_precision_score(Ytest, all_yscores_[key])
            # plot the precision-recall curves
            plt.plot(
                recall,
                precision,
                marker=".",
                label=key + ": " + str(round(auc, 2)) + "_" + str(round(avr_prec, 2)),
            )
            # axis labels
            plt.xlabel("Recall")
            plt.ylabel("Precision")
            # show the legend
            plt.legend()
            # show the plot
    #         plt.show()

    pp = PdfPages(out_fig + "_all.pdf")
    pp.savefig(fig14, bbox_inches="tight")
    pp.savefig(fig15, bbox_inches="tight")
    pp.savefig(fig16, bbox_inches="tight")
    pp.savefig(fig17, bbox_inches="tight")
    pp.savefig(fig18, bbox_inches="tight")
    pp.savefig(fig19, bbox_inches="tight")
    pp.savefig(fig20, bbox_inches="tight")
    pp.savefig(fig21, bbox_inches="tight")

    pp.close()

    # Compute and plot SHAP values
    dim_exp = x_exp.shape[1]
    dim_cells = x_cell.shape[1]

    with open(snakemake.input["training"], "rb") as b:
        training_set = pickle.load(b)
    shap_values_all_exp, shap_values_all_cell, shap_values_all_gender = shap_loop(
        list(model_j), training_set, dim_exp, dim_cells, x_exp, x_cell
    )

    nb = len(model_j)
    f1 = plt.figure()
    with plt.rc_context({"figure.figsize": (4, 3), "figure.dpi": 300}):
        shap.summary_plot(
            shap_values_all_exp / nb,
            plot_type="violin",
            features=np.array(x_exp),
            feature_names=genes,
            color_bar_label="Feature value",
            show=False,
        )
    f2 = plt.figure()

    with plt.rc_context({"figure.figsize": (4, 3), "figure.dpi": 300}):
        shap.summary_plot(
            shap_values_all_exp / nb,
            plot_type="bar",
            features=np.array(x_exp),
            feature_names=genes,
            color_bar_label="Feature value",
            show=False,
        )
    f3 = plt.figure()

    with plt.rc_context({"figure.figsize": (4, 3), "figure.dpi": 300}):
        shap.summary_plot(
            shap_values_all_cell / nb,
            plot_type="violin",
            features=x_cell,
            feature_names=selected_cols,
            color_bar_label="Feature value",
            show=False,
            max_display=15,
        )
    f4 = plt.figure()
    with plt.rc_context({"figure.figsize": (4, 3), "figure.dpi": 300}):
        shap.summary_plot(
            shap_values_all_cell / nb,
            plot_type="bar",
            features=np.array(x_cell),
            feature_names=selected_cols,
            color_bar_label="Feature value",
            show=False,
            max_display=15,
        )

    f5 = plt.figure()
    #     with plt.rc_context({'figure.figsize': (4, 3), 'figure.dpi':300}):
    #         shap.summary_plot(shap_values_all_gender/nb, plot_type= 'bar',features=np.concatenate((gender_dummies,age_), axis=1)
    #                           , feature_names =['Male','Female','Age'],color_bar_label='Feature value',show=False, max_display=15)
    f6 = plt.figure()
    with plt.rc_context({"figure.figsize": (4, 3), "figure.dpi": 300}):
        shap.summary_plot(
            shap_values_all_gender / nb,
            plot_type="bar",
            features=np.concatenate((gender_dummies, age_), axis=1),
            feature_names=["Male", "Female", "Age"],
            color_bar_label="Feature value",
            show=False,
            max_display=15,
        )

    pp = PdfPages(out_shap)
    pp.savefig(f1, bbox_inches="tight")
    pp.savefig(f2, bbox_inches="tight")
    pp.savefig(f3, bbox_inches="tight")
    pp.savefig(f4, bbox_inches="tight")
    pp.savefig(f5, bbox_inches="tight")
    pp.savefig(f6, bbox_inches="tight")
    pp.close()
