"""
Classifiers to differentiate region C and E SNPs (True Positives and False Negatives)
_experiment_= "3"
__author__ = “Jiajia Xu”
__copyright__ = “Copyright 2019, Individual Computing Project”
__email__ = “jiajia.xu@anu.edu.au”
__status__ = “Production”
"""

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.preprocessing import StandardScaler,MinMaxScaler

dataset=pd.read_csv("C_E_dataset.csv",index_col=None)
X=dataset[['qual','repeats','sd','read depth','allele fraction']]
y=dataset['truth']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3,shuffle=True)

# normalization
X_train = pd.DataFrame(MinMaxScaler().fit_transform(X_train), columns = X_train.columns)
X_test = pd.DataFrame(MinMaxScaler().fit_transform(X_test), columns = X_test.columns)

clf=RandomForestClassifier(n_estimators=1000)
clf.fit(X_train,y_train)
y_pred=clf.predict(X_test)

print("Accuracy: ",accuracy_score(y_test,y_pred))

clf1=SVC(gamma='scale')
clf1.fit(X_train,y_train)
y_pred1=clf1.predict(X_test)

print("Accuracy: ",accuracy_score(y_test,y_pred1))

clf2=MLPClassifier(solver='lbfgs', alpha=1e-5,random_state=1)
clf2.fit(X_train,y_train)
y_pred2=clf2.predict(X_test)

print("Accuracy: ",accuracy_score(y_test,y_pred2))

clf3=AdaBoostClassifier(n_estimators=100, random_state=0)
clf3.fit(X_train,y_train)
y_pred3=clf3.predict(X_test)

print("Accuracy: ",accuracy_score(y_test,y_pred3))