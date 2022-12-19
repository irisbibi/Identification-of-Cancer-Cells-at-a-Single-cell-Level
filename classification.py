# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import roc_curve, roc_auc_score, classification_report
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import StratifiedKFold

# %%
with open ('GSE144735_DEGs.txt') as file:
    degs = [line.rstrip() for line in file]
expression_data1 = pd.read_csv('GSE144735_norm_data_rmBorder.csv') 
annotation1 = pd.read_csv('GSE144735_annotation_rmBorder.csv') 
expression_data1_all_features=pd.read_csv('GSE144735_norm_data_rmBorder_allFeat.csv')
expression_data2 = pd.read_csv('GSE200997_norm_data.csv') 
annotation2 = pd.read_csv('GSE200997_annotation.csv')
expression_data3 = pd.read_csv('GSE132465_norm_data_allFeat.csv') 
annotation3 = pd.read_csv('GSE132465_annotation.csv') 

# %%
data1_deg = expression_data1_all_features[expression_data1_all_features.index.isin(degs)].transpose()
deg_expression_with_label = pd.merge(data1_deg, annotation1['Class'], how='left', on=data1_deg.index)
deg_expression_with_label = deg_expression_with_label.set_index('key_0')
deg_expression_with_label.index.names = ['Index']
deg_expression_with_label['Label'] = deg_expression_with_label.Class.replace(to_replace=['Normal', 'Tumor'], value=[0, 1])

# %%
data2_deg = expression_data2[expression_data2.index.isin(degs)].transpose()
deg_expression2 = data2_deg.loc[:,~data2_deg.columns.duplicated()]
deg_expression_with_label2 = pd.merge(deg_expression2, annotation2['Condition'], how='left', on=deg_expression2.index)
deg_expression_with_label2 = deg_expression_with_label2.set_index('key_0')
deg_expression_with_label2.index.names = ['Index']
deg_expression_with_label2['Label'] = deg_expression_with_label2.Condition.replace(to_replace=['Normal', 'Tumor'], value=[0, 1])

# %%
data3_deg = expression_data3[expression_data3.index.isin(degs)].transpose()
deg_expression3 = data3_deg.loc[:,~data3_deg.columns.duplicated()]

deg_expression_with_label3 = pd.merge(deg_expression3, annotation3['Class'], how='left', on=deg_expression3.index)
deg_expression_with_label3 = deg_expression_with_label3.set_index('key_0')
deg_expression_with_label3.index.names = ['Index']
deg_expression_with_label3['Label'] = deg_expression_with_label3.Class.replace(to_replace=['Normal', 'Tumor'], value=[0, 1])

# %%
data1 = deg_expression_with_label.iloc[:,:-2]
label1 = deg_expression_with_label['Label']

# %%
data2 = deg_expression_with_label2.iloc[:,:-2]
label2 = deg_expression_with_label2['Label']

# %%
data3 = deg_expression_with_label3.iloc[:,:-2]
label3 = deg_expression_with_label3['Label']

# %%
x_train, x_test, y_train, y_test = train_test_split(data1, label1, test_size=0.25, random_state=0)
lr1 = LogisticRegression()
lr1.fit(x_train, y_train)
predictions = lr.predict(x_test)
score1 = lr.score(x_test, y_test)
print('Logistic Regression:', score1)
rf1=RandomForestClassifier(n_estimators=30)
rf1.fit(x_train, y_train)
y_pred1_rf=rf1.predict(x_test)
print("Random Forest:", metrics.accuracy_score(y_test, y_pred1_rf))
svm1 = SVC(kernel = 'linear',gamma = 'scale', shrinking = False, probability = True)
svm1.fit(x_train, y_train)
y_pred1_svm = svm1.predict(x_test)
print("SVM: ", metrics.accuracy_score(y_test, y_pred1_svm))

# %%
def roc_curve(X, y, classifier, dataset):
    cv = StratifiedKFold(n_splits=6)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots(figsize=(6, 6))
    for fold, (train, test) in enumerate(cv.split(X, y)):
        classifier.fit(X[train], y[train])
        viz = RocCurveDisplay.from_estimator(
            classifier,
            X[test],
            y[test],
            name=f"ROC fold {fold}",
            alpha=0.3,
            lw=1,
            ax=ax,
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
    ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        xlabel="False Positive Rate",
        ylabel="True Positive Rate",
        title=f"Mean ROC curve of Support Vector Machine for Dataset %d" %dataset,
    )
    ax.axis("square")
    ax.legend(loc="lower right")
    plt.savefig(f'SVM ROC for dataset2.jpg')

# %%
roc_curve(data1.to_numpy(), label1.to_numpy(), lr1,1)

# %%
roc_curve(data1.to_numpy(), label1.to_numpy(), rf1,1)

# %%
roc_curve(data1.to_numpy(), label1.to_numpy(), svm1,1)

# %%
pred_prob1 = lr1.predict_proba(x_test)
pred_prob2 = rf1.predict_proba(x_test)
pred_prob3 = svm1.predict_proba(x_test)
fpr1, tpr1, thresh1 = roc_curve(y_test, pred_prob1[:,1], pos_label=1)
fpr2, tpr2, thresh2 = roc_curve(y_test, pred_prob2[:,1], pos_label=1)
fpr3, tpr3, thresh3 = roc_curve(y_test, pred_prob3[:,1], pos_label=1)
# roc curve for tpr = fpr 
random_probs = [0 for i in range(len(y_test))]
p_fpr, p_tpr, _ = roc_curve(y_test, random_probs, pos_label=1)
lr_auc_score = roc_auc_score(y_test, pred_prob1[:,1])
rf_auc_score = roc_auc_score(y_test, pred_prob2[:,1])
svm_auc_score = roc_auc_score(y_test, pred_prob3[:,1])
print(lr_auc_score)
print(rf_auc_score)
print(svm_auc_score)
plt.style.use('seaborn')

# plot roc curves
plt.plot(fpr1, tpr1, linestyle='--',color='red', label='Logistic Regression')
plt.plot(fpr2, tpr2, linestyle='--',color='orange', label='Random Forest')
plt.plot(fpr3, tpr3, linestyle='--',color='green', label='Support Vector Machine')
plt.plot(p_fpr, p_tpr, linestyle='--', color='blue')
# title
plt.title('ROC curve for Dataset 1')
# x label
plt.xlabel('False Positive Rate')
# y label
plt.ylabel('True Positive rate')

plt.legend(loc='best')
plt.savefig('ROC',dpi=300)
plt.show()

# %%
x_train2, x_test2, y_train2, y_test2 = train_test_split(data2, label2, test_size=0.25, random_state=0)
lr = LogisticRegression()
lr.fit(x_train2, y_train2)
predictions = lr.predict(x_test2)
score = lr.score(x_test2, y_test2)
print('Logistic Regression:', score)
rf=RandomForestClassifier(n_estimators=30)
rf.fit(x_train2, y_train2)
y_pred_data2=rf.predict(x_test2)
print("Random Forest:", metrics.accuracy_score(y_test2, y_pred_data2))
svm = SVC(kernel = 'linear',gamma = 'scale', shrinking = False, probability = True)
svm.fit(x_train2, y_train2)
y_pred_svm_data2 = svm.predict(x_test2)
print("SVM: ", metrics.accuracy_score(y_test2, y_pred_svm_data2))

# %%
roc_curve(data2.to_numpy(), label2.to_numpy(), lr,2)

# %%
roc_curve(data2.to_numpy(), label2.to_numpy(), rf,2)

# %%
roc_curve(data2.to_numpy(), label2.to_numpy(), svm,2)

# %%
pred_prob1 = lr.predict_proba(x_test)
pred_prob2 = rf.predict_proba(x_test)
pred_prob3 = svm.predict_proba(x_test)
fpr1, tpr1, thresh1 = roc_curve(y_test, pred_prob1[:,1], pos_label=1)
fpr2, tpr2, thresh2 = roc_curve(y_test, pred_prob2[:,1], pos_label=1)
fpr3, tpr3, thresh3 = roc_curve(y_test, pred_prob3[:,1], pos_label=1)
# roc curve for tpr = fpr 
random_probs = [0 for i in range(len(y_test))]
p_fpr, p_tpr, _ = roc_curve(y_test, random_probs, pos_label=1)
lr_auc_score = roc_auc_score(y_test, pred_prob1[:,1])
rf_auc_score = roc_auc_score(y_test, pred_prob2[:,1])
svm_auc_score = roc_auc_score(y_test, pred_prob3[:,1])
print(lr_auc_score)
print(rf_auc_score)
print(svm_auc_score)
plt.style.use('seaborn')

# plot roc curves
plt.plot(fpr1, tpr1, linestyle='--',color='red', label='Logistic Regression')
plt.plot(fpr2, tpr2, linestyle='--',color='orange', label='Random Forest')
plt.plot(fpr3, tpr3, linestyle='--',color='green', label='Support Vector Machine')
plt.plot(p_fpr, p_tpr, linestyle='--', color='blue')
# title
plt.title('ROC curve for Dataset 2')
# x label
plt.xlabel('False Positive Rate')
# y label
plt.ylabel('True Positive rate')

plt.legend(loc='best')
plt.savefig('ROC',dpi=300)
plt.show()

# %%
x_train3, x_test3, y_train3, y_test3 = train_test_split(data3, label3, test_size=0.25, random_state=0)
lr3 = LogisticRegression()
lr3.fit(x_train3, y_train3)
predictions = lr.predict(x_test3)
score = lr.score(x_test3, y_test3)
print('Logistic Regression:', score)
rf3 = RandomForestClassifier(n_estimators=30)
rf3.fit(x_train3, y_train3)
y_pred3_rf = rf.predict(x_test3)
print("Random Forest:", metrics.accuracy_score(y_test3, y_pred3_rf))
svm3 = SVC(kernel = 'linear',gamma = 'scale', shrinking = False, probability = True)
svm3.fit(x_train3, y_train3)
y_pred3_svm = svm.predict(x_test3)
print(classification_report(y_test3, y_pred3_svm))

# %%
roc_curve(data3.to_numpy(), label3.to_numpy(), lr3, 3)

# %%
roc_curve(data3.to_numpy(), label3.to_numpy(), rf3, 3)

# %%
roc_curve(data3.to_numpy(), label3.to_numpy(), svm3, 3)

# %%
pred_prob1 = lr.predict_proba(x_test)
pred_prob2 = rf.predict_proba(x_test)
pred_prob3 = svm.predict_proba(x_test)
fpr1, tpr1, thresh1 = roc_curve(y_test, pred_prob1[:,1], pos_label=1)
fpr2, tpr2, thresh2 = roc_curve(y_test, pred_prob2[:,1], pos_label=1)
fpr3, tpr3, thresh3 = roc_curve(y_test, pred_prob3[:,1], pos_label=1)
# roc curve for tpr = fpr 
random_probs = [0 for i in range(len(y_test))]
p_fpr, p_tpr, _ = roc_curve(y_test, random_probs, pos_label=1)
lr_auc_score = roc_auc_score(y_test, pred_prob1[:,1])
rf_auc_score = roc_auc_score(y_test, pred_prob2[:,1])
svm_auc_score = roc_auc_score(y_test, pred_prob3[:,1])
print(lr_auc_score)
print(rf_auc_score)
print(svm_auc_score)
plt.style.use('seaborn')

# plot roc curves
plt.plot(fpr1, tpr1, linestyle='--',color='red', label='Logistic Regression')
plt.plot(fpr2, tpr2, linestyle='--',color='orange', label='Random Forest')
plt.plot(fpr3, tpr3, linestyle='--',color='green', label='Support Vector Machine')
plt.plot(p_fpr, p_tpr, linestyle='--', color='blue')
# title
plt.title('ROC curve for Dataset 3')
# x label
plt.xlabel('False Positive Rate')
# y label
plt.ylabel('True Positive rate')

plt.legend(loc='best')
plt.savefig('ROC',dpi=300)
plt.show()